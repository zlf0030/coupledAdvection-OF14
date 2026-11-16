/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of the levelSetVoF solver module, derived from the
    OpenFOAM-14 incompressibleVoF module (GNU General Public License v3).

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "levelSetTwoPhaseVoFMixture.H"
#include "fvMesh.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "surfaceInterpolate.H"
#include "zeroGradientFvPatchFields.H"
#include "calculatedFvPatchFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(levelSetTwoPhaseVoFMixture, 0);
}


// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

Foam::scalar Foam::levelSetTwoPhaseVoFMixture::minCellSpacing() const
{
    // Cell-centre spacing across internal faces. Unlike cbrt(V) this is
    // also right for 2-D and axisymmetric meshes of any thickness.
    return gMin(mag(levelSet_.mesh().delta())().primitiveField());
}


Foam::levelSetTwoPhaseVoFMixture::couplingSource
Foam::levelSetTwoPhaseVoFMixture::readCouplingSource
(
    const dictionary& coeffs,
    const word& name
) const
{
    const word s(coeffs.lookupOrDefault<word>(name, "levelSet"));

    if (s == "volumeFraction")
    {
        return volumeFraction;
    }

    if (s == "levelSet")
    {
        return levelSet;
    }

    FatalIOErrorInFunction(coeffs)
        << "Unknown " << name << " value " << s
        << "; expected volumeFraction or levelSet"
        << exit(FatalIOError);

    return levelSet;
}


void Foam::levelSetTwoPhaseVoFMixture::readCoupling()
{
    // *this is the phaseProperties dictionary (twoPhaseMixture base)
    const dictionary coeffs(subOrEmptyDict("levelSet"));

    cellSize_ = coeffs.lookupOrDefault<scalar>("cellSize", minCellSpacing());
    widthCoeff_ = coeffs.lookupOrDefault<scalar>("widthCoeff", 1.5);
    seedCoeff_ = coeffs.lookupOrDefault<scalar>("seedCoeff", 0.75);
    pseudoTimeCoeff_ =
        coeffs.lookupOrDefault<scalar>("pseudoTimeCoeff", 0.1);

    if (cellSize_ <= 0 || widthCoeff_ <= 0 || pseudoTimeCoeff_ <= 0)
    {
        FatalIOErrorInFunction(coeffs)
            << "cellSize, widthCoeff and pseudoTimeCoeff must be positive"
            << exit(FatalIOError);
    }

    // Enough pseudo-time to push the distance a full half-width away from
    // the interface, so cells outside the band always see |d| > eps
    nRedistanceIter_ = coeffs.lookupOrDefault<label>
    (
        "nRedistanceIter",
        label(std::ceil(widthCoeff_/pseudoTimeCoeff_ - SMALL))
    );

    curvatureFrom_ = readCouplingSource(coeffs, "curvatureFrom");
    propertiesFrom_ = readCouplingSource(coeffs, "propertiesFrom");
    massImbalance_ = coeffs.lookupOrDefault<bool>("massImbalance", true);

    Info<< "Level-set coupling:" << nl
        << "    cell size       " << cellSize_ << nl
        << "    half-width      " << widthCoeff_*cellSize_ << nl
        << "    redistance      " << nRedistanceIter_ << " x "
        << pseudoTimeCoeff_*cellSize_ << nl
        << "    curvature from  "
        << (curvatureFrom_ == levelSet ? "levelSet" : "volumeFraction") << nl
        << "    properties from "
        << (propertiesFrom_ == levelSet ? "levelSet" : "volumeFraction") << nl
        << "    mass imbalance  "
        << (massImbalance_ ? "on" : "off")
        << nl << endl;
}


Foam::dimensionedScalar Foam::levelSetTwoPhaseVoFMixture::halfWidth() const
{
    return dimensionedScalar(dimLength, widthCoeff_*cellSize_);
}


void Foam::levelSetTwoPhaseVoFMixture::redistance()
{
    const dimensionedScalar h(dimLength, cellSize_);
    const dimensionedScalar dTau(pseudoTimeCoeff_*h);

    const volScalarField seed((2.0*alpha1() - 1.0)*seedCoeff_*h);

    // Sign of the seed; zero on cells sitting exactly on the interface
    const volScalarField seedSign
    (
        seed/max(mag(seed), dimensionedScalar(dimLength, VSMALL))
    );

    levelSet_ == seed;
    levelSet_.correctBoundaryConditions();

    for (label iter = 0; iter < nRedistanceIter_; ++iter)
    {
        levelSet_ += seedSign*(1.0 - mag(fvc::grad(levelSet_)))*dTau;
        levelSet_.correctBoundaryConditions();
    }

    // Report how far |grad d| is from 1 inside the interface band
    const scalarField& V = levelSet_.mesh().V();
    const scalarField inBand
    (
        pos0(halfWidth().value() - mag(levelSet_.primitiveField()))
    );
    const scalarField gradError
    (
        mag(1.0 - mag(fvc::grad(levelSet_)().primitiveField()))
    );
    const scalar bandVolume = gSum(inBand*V);

    if (bandVolume > 0)
    {
        Info<< "Level set: mean |1 - |grad d|| in band = "
            << gSum(inBand*gradError*V)/bandVolume << endl;
    }
}


void Foam::levelSetTwoPhaseVoFMixture::updateHeaviside()
{
    using constant::mathematical::pi;

    // Clamping x to [-1, 1] makes the smooth formula exact outside the band
    const volScalarField x
    (
        max(min(levelSet_/halfWidth(), scalar(1)), scalar(-1))
    );

    // operator== assigns the internal field and the boundary patches
    // alike, so H is evaluated on walls and processor patches (C4)
    heaviside_ == 0.5*(1.0 + x + sin(pi*x)/pi);
}


void Foam::levelSetTwoPhaseVoFMixture::updateCurvature()
{
    const surfaceVectorField gradDf(fvc::interpolate(fvc::grad(levelSet_)));

    const surfaceScalarField nf
    (
        (gradDf/(mag(gradDf) + dimensionedScalar(dimless, 1e-6)))
      & levelSet_.mesh().Sf()
    );

    levelSetCurvature_ = -fvc::div(nf);
}


void Foam::levelSetTwoPhaseVoFMixture::blendProperties()
{
    const volScalarField& H = heaviside_;

    // C1/C3: continuity density from the Heaviside blend. The solver's
    // rho is a const reference bound to the base's private rho_ field at
    // construction, so rewriting the contents of that field here is the
    // only way to make every ddt(rho) and div(rhoPhi) see the blended
    // density without modifying the official solver.
    volScalarField& rho =
        const_cast<volScalarField&>
        (
            incompressibleTwoPhaseVoFMixture::rho()
        );
    rho == H*rho1() + (1.0 - H)*rho2();

    // C3: blend the dynamic viscosity, then divide by the blended
    // density, so that rho*nu is the Heaviside average of mu1 and mu2
    nuBlended_ =
    (
        H*rho1()*nuModel1().nu() + (1.0 - H)*rho2()*nuModel2().nu()
    )/(H*rho1() + (1.0 - H)*rho2());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::levelSetTwoPhaseVoFMixture::rebuildLevelSet()
{
    redistance();
    updateHeaviside();
    updateCurvature();

    if (propertiesFrom_ == levelSet)
    {
        blendProperties();
    }
}


Foam::tmp<Foam::volScalarField> Foam::levelSetTwoPhaseVoFMixture::nu() const
{
    if (propertiesFrom_ == levelSet)
    {
        return nuBlended_;
    }

    return incompressibleTwoPhaseVoFMixture::nu();
}


Foam::tmp<Foam::scalarField> Foam::levelSetTwoPhaseVoFMixture::nu
(
    const label patchi
) const
{
    if (propertiesFrom_ == levelSet)
    {
        return nuBlended_.boundaryField()[patchi];
    }

    return incompressibleTwoPhaseVoFMixture::nu(patchi);
}


void Foam::levelSetTwoPhaseVoFMixture::correct()
{
    incompressibleTwoPhaseVoFMixture::correct();

    if (propertiesFrom_ == levelSet)
    {
        blendProperties();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::levelSetTwoPhaseVoFMixture::levelSetTwoPhaseVoFMixture
(
    const fvMesh& mesh
)
:
    incompressibleTwoPhaseVoFMixture(mesh),
    cellSize_(0),
    widthCoeff_(0),
    seedCoeff_(0),
    pseudoTimeCoeff_(0),
    nRedistanceIter_(0),
    curvatureFrom_(levelSet),
    propertiesFrom_(levelSet),
    massImbalance_(true),
    sigma_("sigma", dimensionSet(1, 0, -2, 0, 0), *this),
    levelSet_
    (
        IOobject
        (
            "levelSet",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimLength, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    heaviside_
    (
        IOobject
        (
            "heaviside",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
    ),
    levelSetCurvature_
    (
        IOobject
        (
            "levelSetCurvature",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    nuBlended_
    (
        IOobject
        (
            "nuBlended",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimKinematicViscosity, Zero),
        calculatedFvPatchScalarField::typeName
    )
{
    readCoupling();
    rebuildLevelSet();
}


// ************************************************************************* //
