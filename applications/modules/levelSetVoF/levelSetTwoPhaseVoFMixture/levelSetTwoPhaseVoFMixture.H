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

Class
    Foam::levelSetTwoPhaseVoFMixture

Description
    Two-phase incompressible VoF mixture with a coupled level-set field
    (port of interLevelSetFoam, Taran, arXiv:2610.24225).

    The level set (signed distance) field is rebuilt from the transported
    volume fraction each time rebuildLevelSet() is called:
      - Seed:             psi0 = (2*alpha1 - 1)*seedCoeff*h        (Eq. 3)
      - Redistance:       dtau*S(psi0)*(1 - |grad psi|) iterations (Eq. 4)
      - Heaviside:        H = 0.5*(1 + x + sin(pi*x)/pi)          (Eq. 5)
      - Curvature:        kappa = -div(nHatf)                     (Eq. 6)

    With curvatureFrom = levelSet the surface tension force uses the
    level-set curvature and the normal gradient of the Heaviside function
    (consistency condition C2).  With propertiesFrom = levelSet the
    continuity density is the Heaviside blend and the kinematic viscosity
    is the Heaviside average of the dynamic viscosities divided by the
    blended density (C3).  Both can be switched back to the volume
    fraction, which gives the volume-fraction coupling.

SourceFiles
    levelSetTwoPhaseVoFMixture.C

\*---------------------------------------------------------------------------*/

#ifndef levelSetTwoPhaseVoFMixture_H
#define levelSetTwoPhaseVoFMixture_H

#include "incompressibleTwoPhaseVoFMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*
              Class levelSetTwoPhaseVoFMixture Declaration
\*---------------------------------------------------------------------------*/

class levelSetTwoPhaseVoFMixture
:
    public incompressibleTwoPhaseVoFMixture
{
public:

    //- Source of the curvature and of the fluid properties
    enum couplingSource
    {
        volumeFraction,
        levelSet
    };


protected:

    // Private Data

        //- Reference cell size h [m]
        scalar cellSize_;

        //- Interface half-width in units of h
        scalar widthCoeff_;

        //- Magnitude of the seed distance in units of h
        scalar seedCoeff_;

        //- Pseudo-time step of the redistance iterations in units of h
        scalar pseudoTimeCoeff_;

        //- Number of redistance iterations
        label nRedistanceIter_;

        //- Source of the curvature in the surface tension force
        couplingSource curvatureFrom_;

        //- Source of the density and viscosity blending
        couplingSource propertiesFrom_;

        //- Include the mass-imbalance correction in the momentum
        //  equation (default true)
        bool massImbalance_;

        //- Surface tension coefficient [kg/s^2] (constant)
        dimensionedScalar sigma_;

        //- Signed distance to the interface [m], positive in phase 1
        volScalarField levelSet_;

        //- Smoothed Heaviside function of the signed distance
        volScalarField heaviside_;

        //- Interface curvature from the signed distance [1/m]
        volScalarField levelSetCurvature_;

        //- Kinematic viscosity blended from the Heaviside function,
        //  used when propertiesFrom = levelSet
        volScalarField nuBlended_;


    // Private Member Functions

        //- Read the coupling settings from the levelSet sub-dictionary
        void readCoupling();

        //- Read a coupling source switch from a dictionary entry
        couplingSource readCouplingSource
        (
            const dictionary& coeffs,
            const word& name
        ) const;

        //- Smallest cell-centre spacing of the mesh
        scalar minCellSpacing() const;

        //- Interface half-width as a dimensioned length
        dimensionedScalar halfWidth() const;

        //- Seed the signed distance from alpha1 and redistance it
        void redistance();

        //- Update the Heaviside function from the signed distance
        void updateHeaviside();

        //- Update the curvature from the signed distance
        void updateCurvature();

        //- Blend the continuity density and the kinematic viscosity
        //  with the Heaviside function (consistency conditions C1, C3)
        void blendProperties();

        //- No copy construct
        levelSetTwoPhaseVoFMixture(const levelSetTwoPhaseVoFMixture&) = delete;

        //- No copy assignment
        void operator=(const levelSetTwoPhaseVoFMixture&) = delete;


public:

    TypeName("levelSetTwoPhaseVoFMixture");


    // Constructors

        //- Construct from a mesh
        levelSetTwoPhaseVoFMixture(const fvMesh& mesh);


    //- Destructor
    virtual ~levelSetTwoPhaseVoFMixture()
    {}


    // Member Functions

        // Access

            //- Signed distance to the interface
            const volScalarField& signedDistance() const
            {
                return levelSet_;
            }

            //- Smoothed Heaviside function
            const volScalarField& heaviside() const
            {
                return heaviside_;
            }

            //- Interface curvature from the signed distance
            const volScalarField& levelSetCurvature() const
            {
                return levelSetCurvature_;
            }

            //- Surface tension coefficient
            const dimensionedScalar& sigma() const
            {
                return sigma_;
            }

            //- Source of the curvature in the surface tension force
            couplingSource curvatureFrom() const
            {
                return curvatureFrom_;
            }

            //- Source of the density and viscosity blending
            couplingSource propertiesFrom() const
            {
                return propertiesFrom_;
            }

            //- Include the mass-imbalance correction?
            bool massImbalance() const
            {
                return massImbalance_;
            }


        // Coupling

            //- Rebuild the level set, Heaviside, curvature and (if
            //  propertiesFrom = levelSet) the blended properties from
            //  the current alpha1
            void rebuildLevelSet();

        // Viscosity (overrides)

            //- Mixture kinematic viscosity
            virtual tmp<volScalarField> nu() const;

            //- Mixture kinematic viscosity for a patch
            virtual tmp<scalarField> nu(const label patchi) const;

            //- Correct the mixture properties
            virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
