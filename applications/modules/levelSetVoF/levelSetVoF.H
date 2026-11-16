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
    Foam::solvers::levelSetVoF

Description
    Solver module for 2 incompressible, isothermal immiscible fluids using a
    coupled level-set/VOF (CLSVOF) interface capturing approach, with optional
    mesh motion and mesh topology changes including adaptive re-meshing.

    Stage-1 shell: structurally identical to the incompressibleVoF module
    except that the mixture is a levelSetTwoPhaseVoFMixture.  Stage 2 added
    the level-set field reconstruction (seed -> redistance -> Heaviside ->
    curvature).  Stage 3 couples it into the momentum equation and the
    material properties:
      - C2: surface tension force from sigma*kappa*snGrad(H)
      - C1: mass-imbalance correction in the momentum predictor
      - C3: density and viscosity blended from the Heaviside function

    Reference: Taran, A. (2026), "A consistent coupled level-set and
    volume-of-fluid solver for two-phase flow in OpenFOAM", arXiv:2610.24225

SourceFiles
    levelSetVoF.C

See also
    Foam::solvers::VoFSolver
    Foam::solvers::twoPhaseVoFSolver
    Foam::solvers::incompressibleVoF

\*---------------------------------------------------------------------------*/

#ifndef levelSetVoF_H
#define levelSetVoF_H

#include "twoPhaseVoFSolver.H"
#include "levelSetTwoPhaseVoFMixture.H"
#include "incompressibleInterPhaseTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*
                           Class levelSetVoF Declaration
\*---------------------------------------------------------------------------*/

class levelSetVoF
:
    public twoPhaseVoFSolver
{

protected:

    // Phase properties

        //- The incompressible two-phase mixture with level-set support
        levelSetTwoPhaseVoFMixture& mixture;


    // Thermophysical properties

        //- Static pressure field
        volScalarField p;


    // Pressure reference

        //- Pressure reference
        Foam::pressureReference pressureReference_;


    // Momentum transport

        //- Momentum transport model
        incompressibleInterPhaseTransportModel momentumTransport;


    // Protected Member Functions

        //- Return the pressure reference
        virtual const Foam::pressureReference& pressureReference() const
        {
            return pressureReference_;
        }

        //- The flow is incompressible
        virtual bool incompressible() const
        {
            return true;
        }

        //- Is the flow divergent?
        //  i.e. includes phase-fraction sources
        virtual bool divergent() const
        {
            return
                fvModels().addsSupToField(alpha1.name())
             || fvModels().addsSupToField(alpha2.name());
        }

        //- Return the mixture compressibility/density
        //  Not required for incompressible fluids
        virtual tmp<volScalarField> psiByRho() const
        {
            return tmp<volScalarField>(nullptr);
        }

        //- Calculate the alpha equation sources
        virtual void alphaSuSp
        (
            tmp<volScalarField::Internal>& Su,
            tmp<volScalarField::Internal>& Sp
        );

        //- Return the interface surface tension force for the momentum
        //  equation; uses the level-set curvature and the normal gradient
        //  of the Heaviside function when curvatureFrom = levelSet (C2)
        virtual tmp<surfaceScalarField> surfaceTensionForce() const;

        //- Momentum predictor with the mass-imbalance correction that
        //  keeps the momentum equation consistent with the blended
        //  density and mass flux (C1)
        virtual void momentumPredictor();

        //- Return the momentum equation stress term
        virtual tmp<fvVectorMatrix> divDevTau(volVectorField& U)
        {
            return momentumTransport.divDevTau(rho, U);
        }


public:

    //- Runtime type information
    TypeName("levelSetVoF");


    // Constructors

        //- Construct from region mesh
        levelSetVoF(fvMesh& mesh);

        //- Construct from mesh and a custom two-phase VoF mixture
        //  (must derive from levelSetTwoPhaseVoFMixture)
        levelSetVoF
        (
            fvMesh& mesh,
            autoPtr<twoPhaseVoFMixture> mixturePtr
        );

        //- Disallow default bitwise copy construction
        levelSetVoF(const levelSetVoF&) = delete;


    //- Destructor
    virtual ~levelSetVoF();


    // Member Functions

        //- Called at the start of the PIMPLE loop
        virtual void prePredictor();

        //- Predict the momentum transport
        virtual void momentumTransportPredictor();

        //- Predict thermophysical transport
        virtual void thermophysicalTransportPredictor();

        //- Construct and solve the energy equation,
        //  convert to temperature
        //  and update thermophysical and transport properties
        virtual void thermophysicalPredictor();

        //- Construct and solve the pressure equation in the PISO loop
        virtual void pressureCorrector();

        //- Correct the momentum transport
        virtual void momentumTransportCorrector();

        //- Correct the thermophysical transport
        virtual void thermophysicalTransportCorrector();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const levelSetVoF&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
