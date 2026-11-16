/*---------------------------------------------------------------------------*\
  plicLevelSetVoF source
\*---------------------------------------------------------------------------*/

#include "plicLevelSetVoF.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceInterpolate.H"
#include "fvcGrad.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(plicLevelSetVoF, 0);
    addToRunTimeSelectionTable(solver, plicLevelSetVoF, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::plicLevelSetVoF::plicLevelSetVoF(fvMesh& mesh)
:
    levelSetVoF
    (
        mesh,
        autoPtr<twoPhaseVoFMixture>
        (
            new plicLevelSetTwoPhaseVoFMixture(mesh)
        )
    ),
    advectNormalFromLevelSet_
    (
        mesh.solution().solverDict(alpha1.name()).lookupOrDefault<word>
        (
            "advectNormal",
            "default"
        ) == "fromLevelSet"
    ),
    nHatLS_
    (
        IOobject
        (
            "nHatLS",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector(dimless, Zero),
        zeroGradientFvPatchVectorField::typeName
    ),
    plicMixture
    (
        refCast<plicLevelSetTwoPhaseVoFMixture>(twoPhaseVoFSolver::mixture)
    ),
    advector(new isoAdvection(alpha1, phi, U))
{
    // Redo the level-set rebuild with the PLIC-geometry seed (the mixture
    // constructor used the alpha-approximation seed before the flag was
    // available)
    if (plicMixture.seedPlic())
    {
        plicMixture.rebuildLevelSet();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::plicLevelSetVoF::~plicLevelSetVoF()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::plicLevelSetVoF::prePredictor()
{
    // PLIC cutting-plane orientation. When enabled, use the level set
    // rebuilt during the previous step (semantics of the original
    // coupled LS/VOF chain: advect with the previous LS normal).
    if (advectNormalFromLevelSet_)
    {
        const volVectorField gradPsi(fvc::grad(mixture.signedDistance()));
        nHatLS_ =
            gradPsi/max(mag(gradPsi), dimensionedScalar(dimless, SMALL));

        advector->setExternalNormal(nHatLS_);
    }
    else
    {
        advector->clearExternalNormal();
    }

    // PLIC alpha transport over the whole time-step; isoAdvection handles
    // its own sub-cycling. Replaces twoPhaseSolver::alphaPredictor (MULES).
    advector->advect();

    // Alpha fluxes from the face volume transport dVf [m^3]
    alphaPhi1 = advector->getDvf()/runTime.deltaT();
    alpha2 = scalar(1) - alpha1;
    alphaPhi2 = phi - alphaPhi1;

    // Refresh the interface geometry cache and mixture properties, as
    // twoPhaseSolver::prePredictor (alphaPredictor + mixture.correct) would
    correctInterface();
    mixture.correct();

    // Rebuild the level set from the freshly transported alpha, then
    // apply the C1 mass flux handling of levelSetVoF::prePredictor
    plicMixture.rebuildLevelSet();

    if
    (
        mixture.propertiesFrom()
     == levelSetTwoPhaseVoFMixture::couplingSource::levelSet
    )
    {
        // C1: mass flux consistent with the Heaviside-blended density
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

        rhoPhi = alphaPhi1*rho1 + alphaPhi2*rho2;
    }
}


// ************************************************************************* //
