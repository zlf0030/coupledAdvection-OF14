/*---------------------------------------------------------------------------*\
  plicLevelSetVoF - PLIC alpha transport + rebuilt level-set coupling

  Stage D: inherits levelSetVoF (redistance/Heaviside/curvature chain,
  C1-C4 consistency conditions, sigma*kappa*snGrad(H) surface tension)
  and replaces the MULES alpha transport of
  twoPhaseSolver::alphaPredictor with isoAdvection PLIC advection in
  prePredictor().

  Stage E: optional LS-driven PLIC cutting normals (advectNormal
  fromLevelSet) using the external-normal mode of isoAdvection; the
  normal is taken from the level set rebuilt during the previous step.
\*---------------------------------------------------------------------------*/

#ifndef plicLevelSetVoF_H
#define plicLevelSetVoF_H

#include "levelSetVoF.H"
#include "isoAdvection.H"
#include "plicLevelSetTwoPhaseVoFMixture.H"

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                        Class plicLevelSetVoF Declaration
\*---------------------------------------------------------------------------*/

class plicLevelSetVoF
:
    public levelSetVoF
{

protected:

    // Protected data

        //- PLIC cutting-plane orientation from the previous-step level
        //  set? (solverDict entry: advectNormal default|fromLevelSet)
        bool advectNormalFromLevelSet_;

        //- Interface normal field derived from the level set
        volVectorField nHatLS_;

        //- The mixture with PLIC-geometry level-set seeding
        plicLevelSetTwoPhaseVoFMixture& plicMixture;

        //- PLIC advector operating on alpha1/phi/U (references)
        autoPtr<isoAdvection> advector;


public:

    //- Runtime type information
    TypeName("plicLevelSetVoF");

    // Constructors

        //- Construct from mesh
        plicLevelSetVoF(fvMesh& mesh);

        //- Disallow copy construct
        plicLevelSetVoF(const plicLevelSetVoF&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const plicLevelSetVoF&) = delete;


    //- Destructor
    virtual ~plicLevelSetVoF();


    // Member Functions

        //- Pre-predictor: optional LS normal setup, PLIC alpha transport
        //  (replaces the MULES alphaPredictor), then the inherited
        //  level-set rebuild and the C1 rhoPhi handling of
        //  levelSetVoF::prePredictor
        virtual void prePredictor();
};

} // End namespace solvers
} // End namespace Foam

#endif

// ************************************************************************* //
