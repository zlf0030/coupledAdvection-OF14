/*---------------------------------------------------------------------------*\
  plicLevelSetTwoPhaseVoFMixture - level-set mixture with PLIC-geometry
  seeding of the signed distance (Stage E).

  The signed-distance seed in band cells and their 1-ring neighbourhood
  is taken as the signed distance to the nearest PLIC interface plane
  (extracted with isoCutCell from the current alpha field) instead of
  the alpha approximation (2*alpha1 - 1)*seedCoeff*h. Cells beyond the
  1-ring keep the alpha-approximation seed. Dict switch (phaseProperties
  -> levelSet sub-dictionary): seedPlic (default no).
\*---------------------------------------------------------------------------*/

#ifndef plicLevelSetTwoPhaseVoFMixture_H
#define plicLevelSetTwoPhaseVoFMixture_H

#include "levelSetTwoPhaseVoFMixture.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class plicLevelSetTwoPhaseVoFMixture Declaration
\*---------------------------------------------------------------------------*/

class plicLevelSetTwoPhaseVoFMixture
:
    public levelSetTwoPhaseVoFMixture
{
protected:

    // Protected data

        //- Seed the signed distance from the PLIC interface geometry?
        bool seedPlic_;


public:

    //- Runtime type information
    TypeName("plicLevelSetTwoPhaseVoFMixture");

    // Constructors

        //- Construct from mesh
        explicit plicLevelSetTwoPhaseVoFMixture(const fvMesh& mesh);

        //- Disallow copy construct
        plicLevelSetTwoPhaseVoFMixture
        (
            const plicLevelSetTwoPhaseVoFMixture&
        ) = delete;


    //- Destructor
    virtual ~plicLevelSetTwoPhaseVoFMixture()
    {}


    // Member Functions

        //- PLIC seeding enabled?
        bool seedPlic() const
        {
            return seedPlic_;
        }

        //- Rebuild the level set (PLIC-geometry seed if enabled),
        //  Heaviside, curvature and (if propertiesFrom = levelSet) the
        //  blended properties from the current alpha1
        void rebuildLevelSet();
};

} // End namespace Foam

#endif

// ************************************************************************* //
