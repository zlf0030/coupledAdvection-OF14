/*---------------------------------------------------------------------------*\
  coupledAdvection - PLIC-VoF solver module for OpenFOAM v14
  Stage C: prePredictor() replaces the MULES alpha transport of
  twoPhaseSolver::alphaPredictor with isoAdvection-based PLIC advection.
\*---------------------------------------------------------------------------*/

#ifndef coupledAdvection_H
#define coupledAdvection_H

#include "incompressibleVoF.H"
#include "isoAdvection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                        Class coupledAdvection Declaration
\*---------------------------------------------------------------------------*/

class coupledAdvection
:
    public incompressibleVoF
{

protected:

    //- PLIC advector operating on alpha1/phi/U (references)
    autoPtr<isoAdvection> advector;


public:

    //- Runtime type information
    TypeName("coupledAdvection");

    // Constructors

        //- Construct from mesh
        coupledAdvection(fvMesh& mesh);

        //- Disallow copy construct
        coupledAdvection(const coupledAdvection&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const coupledAdvection&) = delete;


    //- Destructor
    virtual ~coupledAdvection();


    // Member Functions

        //- Pre-predictor: PLIC alpha transport.
        //  Replaces twoPhaseSolver::alphaPredictor (MULES) with
        //  isoAdvection::advect(); derives alphaPhi1 from the face volume
        //  transport dVf, then rhoPhi, then refreshes mixture properties
        //  and the interface geometry caches.
        virtual void prePredictor();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

#endif

// ************************************************************************* //
