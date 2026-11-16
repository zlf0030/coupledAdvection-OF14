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

#include "levelSetVoF.H"
#include "localEulerDdtScheme.H"
#include "fvCorrectPhi.H"
#include "geometricZeroField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(levelSetVoF, 0);
    addToRunTimeSelectionTable(solver, levelSetVoF, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::levelSetVoF::levelSetVoF(fvMesh& mesh)
:
    levelSetVoF
    (
        mesh,
        autoPtr<twoPhaseVoFMixture>(new levelSetTwoPhaseVoFMixture(mesh))
    )
{}


Foam::solvers::levelSetVoF::levelSetVoF
(
    fvMesh& mesh,
    autoPtr<twoPhaseVoFMixture> mixturePtr
)
:
    twoPhaseVoFSolver(mesh, mixturePtr),

    mixture
    (
        refCast<levelSetTwoPhaseVoFMixture>(twoPhaseVoFSolver::mixture)
    ),

    p
    (
        IOobject
        (
            "p",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rho*buoyancy.gh
    ),

    pressureReference_
    (
        p,
        p_rgh,
        pimple.dict()
    ),

    momentumTransport
    (
        U,
        phi,
        alphaPhi1,
        alphaPhi2,
        mixture
    )
{
    if (correctPhi || mesh.topoChanging())
    {
        rAU = new volScalarField
        (
            IOobject
            (
                "rAU",
                runTime.name(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimTime/dimDensity, 1)
        );
    }

    if (!runTime.restart() || !divergent())
    {
        correctUphiBCs(U_, phi_, true);

        fv::correctPhi
        (
            phi_,
            U,
            p_rgh,
            rAU,
            autoPtr<volScalarField>(),
            pressureReference(),
            pimple
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::levelSetVoF::~levelSetVoF()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::levelSetVoF::prePredictor()
{
    twoPhaseVoFSolver::prePredictor();

    // The level set is rebuilt from the freshly transported alpha after
    // the mixture properties have been corrected (cf. interLevelSetFoam:
    // mixture.correct(); mixture.rebuildLevelSet(); before UEqn).
    mixture.rebuildLevelSet();

    if
    (
        mixture.propertiesFrom()
     == levelSetTwoPhaseVoFMixture::couplingSource::levelSet
    )
    {
        // C1: mass flux consistent with the Heaviside-blended density
        // (cf. blendDensity.H in the reference)
        rhoPhi =
            phi_
           *(
                mixture.rho2()
              + (mixture.rho1() - mixture.rho2())
               *fvc::interpolate(mixture.heaviside())
            );
    }
    else
    {
        const dimensionedScalar& rho1 = mixture.rho1();
        const dimensionedScalar& rho2 = mixture.rho2();

        // Calculate the mass-flux
        rhoPhi = alphaPhi1*rho1 + alphaPhi2*rho2;
    }
}


Foam::tmp<Foam::surfaceScalarField>
Foam::solvers::levelSetVoF::surfaceTensionForce() const
{
    if
    (
        mixture.curvatureFrom()
     == levelSetTwoPhaseVoFMixture::couplingSource::levelSet
    )
    {
        // C2: interpolated sigma*kappa times the normal gradient of the
        // Heaviside function (cf. capillaryForce in the reference). Used
        // identically in the momentum predictor and the pressure equation.
        return
            fvc::interpolate
            (
                mixture.sigma()*mixture.levelSetCurvature()
            )*fvc::snGrad(mixture.heaviside());
    }

    return twoPhaseVoFSolver::surfaceTensionForce();
}


void Foam::solvers::levelSetVoF::momentumPredictor()
{
    volVectorField& U = U_;

    // C1: mass imbalance of the rho and rhoPhi used below. Taking it out
    // keeps the momentum equation consistent when rho comes from the
    // level set and is not transported by rhoPhi exactly; it is
    // negligible for the volume-fraction coupling. Can be switched off
    // with the massImbalance entry of the levelSet dictionary.
    if (mixture.massImbalance())
    {
        const volScalarField massImbalance(fvc::ddt(rho) + fvc::div(rhoPhi));

        tUEqn =
        (
            fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
          + fvm::SuSp(-massImbalance, U)
          + MRF.DDt(rho, U)
          + divDevTau(U)
         ==
            fvModels().source(rho, U)
        );
    }
    else
    {
        tUEqn =
        (
            fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
          + MRF.DDt(rho, U)
          + divDevTau(U)
         ==
            fvModels().source(rho, U)
        );
    }
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvConstraints().constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    surfaceTensionForce()
                  - buoyancy.ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );

        fvConstraints().constrain(U);
    }
}


void Foam::solvers::levelSetVoF::momentumTransportPredictor()
{
    momentumTransport.predict();
}


void Foam::solvers::levelSetVoF::thermophysicalTransportPredictor()
{}


void Foam::solvers::levelSetVoF::pressureCorrector()
{
    incompressiblePressureCorrector(p);
}


void Foam::solvers::levelSetVoF::thermophysicalPredictor()
{}


void Foam::solvers::levelSetVoF::momentumTransportCorrector()
{
    momentumTransport.correct();
}


void Foam::solvers::levelSetVoF::thermophysicalTransportCorrector()
{}


// ************************************************************************* //
