/*---------------------------------------------------------------------------*\
  plicLevelSetTwoPhaseVoFMixture source
\*---------------------------------------------------------------------------*/

#include "plicLevelSetTwoPhaseVoFMixture.H"
#include "isoCutCell.H"
#include "volPointInterpolation.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(plicLevelSetTwoPhaseVoFMixture, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::plicLevelSetTwoPhaseVoFMixture::plicLevelSetTwoPhaseVoFMixture
(
    const fvMesh& mesh
)
:
    levelSetTwoPhaseVoFMixture(mesh),
    seedPlic_(subOrEmptyDict("levelSet").lookupOrDefault<bool>("seedPlic", false))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::plicLevelSetTwoPhaseVoFMixture::rebuildLevelSet()
{
    if (!seedPlic_)
    {
        levelSetTwoPhaseVoFMixture::rebuildLevelSet();
        return;
    }

    const fvMesh& mesh = levelSet_.mesh();
    const dimensionedScalar h(dimLength, cellSize_);
    const dimensionedScalar dTau(pseudoTimeCoeff_*h);

    // --- PLIC interface planes from the current alpha field ---

    pointScalarField ap(volPointInterpolation::New(mesh).interpolate(alpha1()));

    DynamicList<label> planeCells;
    DynamicList<vector> planeNormals;    // unit, pointing to the alpha1-low side
    DynamicList<point>  planeCentres;
    {
        isoCutCell cutCell(mesh, ap);

        const volScalarField& a = alpha1();
        forAll(a, celli)
        {
            const label status =
                cutCell.vofCutCell(celli, a[celli], 1e-8, 100);

            if (status == 0)
            {
                const vector n(cutCell.isoFaceArea());
                if (mag(n) > ROOTVSMALL)
                {
                    planeCells.append(celli);
                    planeNormals.append(n/mag(n));
                    planeCentres.append(cutCell.isoFaceCentre());
                }
            }
        }
    }

    const label nPlanes = planeCells.size();

    // Map: cell -> index of its own PLIC plane (-1 if not a surface cell)
    labelList planeOfCell(mesh.nCells(), -1);
    forAll(planeCells, i)
    {
        planeOfCell[planeCells[i]] = i;
    }

    // --- Seed: signed distance to the nearest PLIC plane for surface
    //           cells and their 1-ring neighbourhood, clipped to the
    //           alpha-seed saturation level; alpha approximation
    //           elsewhere ---

    volScalarField seed
    (
        IOobject
        (
            "plicSeed",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (2.0*alpha1() - 1.0)*seedCoeff_*h
    );

    if (nPlanes)
    {
        const vectorField& centres = mesh.cellCentres();
        const labelListList& cellCells = mesh.cellCells();

        forAll(seed, celli)
        {
            scalar bestDist = GREAT;
            scalar bestSigned = 0;

            // Candidate planes: the cell itself and its 1-ring
            if (planeOfCell[celli] != -1)
            {
                const label i = planeOfCell[celli];
                const scalar d =
                    -((centres[celli] - planeCentres[i]) & planeNormals[i]);
                if (mag(d) < bestDist)
                {
                    bestDist = mag(d);
                    bestSigned = d;
                }
            }

            forAll(cellCells[celli], ni)
            {
                const label nei = cellCells[celli][ni];
                if (planeOfCell[nei] != -1)
                {
                    const label i = planeOfCell[nei];
                    const scalar d =
                        -((centres[celli] - planeCentres[i]) & planeNormals[i]);
                    if (mag(d) < bestDist)
                    {
                        bestDist = mag(d);
                        bestSigned = d;
                    }
                }
            }

            if (bestDist < GREAT)
            {
                // Accept the geometric distance only when its sign agrees
                // with the alpha approximation. At high-curvature regions a
                // neighbouring cell can fall on the wrong side of a PLIC
                // plane; a few such sign-flipped seeds each produce an
                // O(1/h) curvature spike, so they are replaced by the
                // always-consistent alpha seed. The clip keeps distances
                // beyond seedCoeff*h at the alpha-seed plateau level so
                // that all seed gradients stay below unity.
                const scalar aSeed = (2.0*alpha1()[celli] - 1.0)*seedCoeff_*cellSize_;
                if
                (
                    sign(bestSigned) == sign(aSeed)
                 || mag(aSeed) < 0.5*seedCoeff_*cellSize_
                )
                {
                    seed[celli] =
                        sign(bestSigned)
                       *min(mag(bestSigned), seedCoeff_*cellSize_);
                }
            }
        }
    }

    // --- Pseudo-time redistance (as levelSetTwoPhaseVoFMixture::redistance)

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

    updateHeaviside();
    updateCurvature();

    if (propertiesFrom() == couplingSource::levelSet)
    {
        blendProperties();
    }
}


// ************************************************************************* //
