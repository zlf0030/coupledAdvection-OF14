/*---------------------------------------------------------------------------*\
  coupledAdvection - PLIC-VoF solver module for OpenFOAM v14
\*---------------------------------------------------------------------------*/

#include "coupledAdvection.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(coupledAdvection, 0);
    addToRunTimeSelectionTable(solver, coupledAdvection, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::coupledAdvection::coupledAdvection(fvMesh& mesh)
:
    incompressibleVoF(mesh),
    advector
    (
        new isoAdvection(alpha1, phi, U)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::coupledAdvection::~coupledAdvection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solvers::coupledAdvection::prePredictor()
{
    // PLIC alpha transport: advance alpha1 over the whole time-step and
    // compute the face volume transport dVf [m^3]. Replaces the MULES
    // alphaSolve of twoPhaseSolver::prePredictor.
    advector->advect();

    // Alpha fluxes from the face volume transport: alphaPhi = dVf/deltaT
    alphaPhi1 = advector->getDvf()/runTime.deltaT();
    alpha2 = scalar(1) - alpha1;
    alphaPhi2 = phi - alphaPhi1;

    // Density flux, as in incompressibleVoF::prePredictor
    rhoPhi = alphaPhi1*mixture.rho1() + alphaPhi2*mixture.rho2();

    // Refresh mixture properties and interface geometry caches
    // (nHatf_/K_ are derived from the updated alpha1 by interface.correct())
    mixture.correct();
    correctInterface();
}


// ************************************************************************* //
