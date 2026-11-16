# coupledAdvection-OF14

A suite of two-phase incompressible interface-capturing solver modules for
**OpenFOAM v14** (OpenFOAM Foundation modular-solver architecture), combining
geometric PLIC advection (isoAdvector) with a rebuilt level-set (LS) field for
interface curvature and physical properties.

Four solvers cover the classic 2x2 combinations:

| Solver (`controlDict` `solver` key) | alpha transport | Curvature / properties |
|---|---|---|
| `incompressibleVoF` (stock OF14) | MULES | alpha-gradient |
| `coupledAdvection` | PLIC (isoAdvection) | alpha-gradient |
| `levelSetVoF` | MULES | rebuilt level set |
| `plicLevelSetVoF` | PLIC (isoAdvection) | rebuilt level set |

The flagship module `plicLevelSetVoF` couples the geometric PLIC transport with
a redistanced level-set field (Heaviside-based curvature and properties) and
implements the consistency conditions of the CLSVOF literature (mass/volume
consistency, zero-level-set fidelity).

## Repository layout

```
src/isoAdvection/                 isoAdvection library ported to OpenFOAM v14
applications/modules/
    coupledAdvection/             modular solver: pure PLIC (VOF-PLIC)
    levelSetVoF/                  modular solver: MULES + level set (CLSVOF-MULES)
    plicLevelSetVoF/              modular solver: PLIC + level set (CLSVOF-PLIC)
cases/hysingBench/                Hysing (2009) rising-bubble benchmark case template
                                  + metric extraction tools + run/clean scripts
Allwmake                          build everything in dependency order
```

## Prerequisites

- OpenFOAM v14 (OpenFOAM Foundation, https://openfoam.org), environment sourced.
- A C++ toolchain supported by OpenFOAM (gcc >= 7).

## Build

The `Make/options` files locate the in-repo libraries through
`$WM_PROJECT_USER_DIR`, so point it at the repository root before building:

```bash
source /path/to/OpenFOAM-v14/etc/bashrc

export WM_PROJECT_USER_DIR=$PWD                 # repo root
export FOAM_USER_LIBBIN=$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib

./Allwmake
```

Build order (handled by `./Allwmake`):
`src/isoAdvection` -> `coupledAdvection` -> `levelSetVoF` -> `plicLevelSetVoF`.

Alternatively, place (or symlink) `src/` and `applications/` into your existing
`$WM_PROJECT_USER_DIR` and build each module with `wmake libso`.

## Run: Hysing rising-bubble benchmark

```bash
cd cases/hysingBench
# 1. generate the exact-Gauss volume-fraction initial field (optional; a
#    pre-generated alpha.air at 1/h=160 is included in 0/)
python3 gen_alpha.py 160
# 2. pick a solver in system/controlDict, e.g.
#        solver plicLevelSetVoF;
# 3. run (4-proc parallel with logging, then reconstruct)
./runPar.sh            # or: nohup ./runPar.sh > log.runPar 2>&1 &
# 4. extract benchmark metrics
python3 metrics.py $PWD 160 320
```

`runPar.sh` / `cleanPar.sh` assume they are executed inside the case directory.

### Solver configuration switches (`plicLevelSetVoF`)

In `system/fvSolution`, solver sub-dictionary `"alpha.air.*"`:

- `advectNormal default|fromLevelSet` - normal used by the PLIC interface
  reconstruction during advection: `default` uses the point-interpolated
  alpha gradient (isoAdvector original); `fromLevelSet` uses the gradient of
  the rebuilt level-set field (strong coupling). Both modes are stable in
  h80/h160 when the time step is resolution-appropriate (deltaT = 0.25*h/U).
- In `constant/phaseProperties`, `levelSet` sub-dictionary:
  `seedPlic on|off` (opt-in, default `off`) - seeds the SDF redistancing from
  PLIC interface geometry instead of the alpha approximation.

## Validation results

Hysing et al. (2009) test 2 (rising bubble), 1/h = 160 (160x320 mesh),
fixed deltaT = 3.90625e-04, endTime = 3 s, 4 MPI ranks.
Reference row: Taran et al., Table 2 (volume-fraction coupling, 1/h = 160).

| Solver | Vc,1 | Vc,2 | yc(3s) | bubble volume drift |
|---|---|---|---|---|
| Taran et al. (paper) | 0.2472 | 0.2315 | 1.1174 | -0.0003% |
| VOF-MULES (`incompressibleVoF`) | 0.2471 | 0.2309 | 1.1142 | -0.0003% |
| CLSVOF-MULES (`levelSetVoF`) | 0.2472 | 0.2315 | 1.1174 | -0.0003% |
| VOF-PLIC (`coupledAdvection`) | 0.2466 | 0.2303 | 1.1117 | -0.0000% |
| CLSVOF-PLIC (`plicLevelSetVoF`) | 0.2471 | 0.2310 | 1.1165 | -0.0000% |
| CLSVOF-PLIC strong coupling (`advectNormal fromLevelSet`) | 0.2471 | 0.2309 | 1.1162 | -0.0000% |

Static droplet parasitic-current test (80x80, t = 1 s, max|U|):

| Curvature source | max|U| |
|---|---|
| alpha-gradient (pure PLIC) | 4.4e-01 |
| level set (MULES + LS) | 9.9e-02 |
| level set (PLIC + LS) | 1.2e-01 |

The level-set curvature reduces spurious currents by roughly 3.8x compared to
the alpha-gradient formulation.

## Provenance and licensing

This entire work is a derivative of GPL-licensed software and is therefore
distributed under the **GNU General Public License v3** (see `LICENSE`).

- `src/isoAdvection` is a port of the
  [isoAdvector](https://github.com/isoAdvector/isoAdvector) library
  (Roenby, Bredmose, Jasak, "A computational method for sharp interface
  advection", J. Comput. Phys. 2016) to the OpenFOAM v14 modular solver
  architecture. Original copyright notices are preserved in the source files;
  the v14 port is marked as modified work.
- All solver modules link against OpenFOAM v14 (GPL-3.0, OpenFOAM Foundation).
- `plicLevelSetVoF` extends `levelSetVoF`, which implements the rebuilt
  level-set coupling chain of the [coupledAdvection](https://github.com/zlf0030/coupledAdvection)
  project, restructured for OpenFOAM v14.

## References

1. Roenby J, Bredmose H, Jasak H. A computational method for sharp interface
   advection. Royal Society Open Science, 2016.
2. Hysing S, Turek S, Kuzmin D, et al. Quantitative benchmark computations of
   two-dimensional bubble dynamics. IJNMF, 2009.
3. Taran F, et al. CLSVOF solver benchmark data (Table 2, volume-fraction
   coupling).
4. OpenFOAM v14 (OpenFOAM Foundation) documentation - modular solvers.
   https://openfoam.org
