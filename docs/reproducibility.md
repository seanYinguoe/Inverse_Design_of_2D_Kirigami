# Running and extending the code

## Requirements

| Task | Requirements |
|---|---|
| Deployment demo and geometry checks | MATLAB; checked with R2025a |
| Inverse design | MATLAB + Optimization Toolbox (`fmincon`) |
| Optional parallel solver setting | Parallel Computing Toolbox |
| Mechanical refinement | COMSOL Multiphysics, appropriate structural-mechanics capability and LiveLink for MATLAB |
| Some legacy fitting / symbolic scripts | Curve Fitting Toolbox / Symbolic Math Toolbox |

Minimum supported MATLAB and COMSOL versions have not been established. The new examples were checked on R2025a; no Octave compatibility is claimed.

## A repeatable run

1. Open this repository root in MATLAB and run `setup_project`.
2. Run `demo_inverse_design('rigid')` or `demo_inverse_design('nonrigid')`.
3. Open the new `results/.../inverse_design.mat`. It contains configuration, solver options, geometry, MATLAB version, runtime and feasibility diagnostics.
4. Retain the Git commit with the run when using it in a manuscript. The full solver reference is in `kinematics-reference.md`; it predates this reorganisation, so its old folder names should be read using the mapping below.

Both grid dimensions should be even for the inherited symmetry machinery. Example dimensions are 4 × 4 units, each containing four panels. Example coordinates use a consistent arbitrary length scale; they are not automatically millimetres.

## Mechanical refinement

The original `GGA*.m` scripts are preserved in `legacy/comsol/`. They are not interchangeable: some optimise full node lists, others exploit a symmetric quarter of the pattern. They assume a `tessellation_compacted` variable and contain their own dimensions, target definitions, material constants and applied displacement. Target selectors also differ from the stage-1 solver. Some final plotting blocks reference additional workspace variables such as `x_initial`.

To work on these drivers, start MATLAB connected to your licensed COMSOL server, run `setup_project`, and add `src/comsol` to the path. Select one case, supply its matching compact geometry, and review the driver and corresponding `livelink*.m` together. Do not substitute the quick-start result without reconciling grid dimensions and boundary conditions.

The cleanup removes personal machine paths from model output locations. It does not change the inherited physical constants or validate every COMSOL model. Long GA/FE runs and fabrication outputs remain outside the verified quick start.

## Exporting geometry

`create_svg_tessellation` writes to `results/exports/`. Its default output is a slit pattern; a nonzero cut width requests a finite-width void. Its hinge detection assumes a regular, compact rotating-squares topology. The quick-start does not automatically export the nonuniform optimised result as a fabrication-ready file.

## Where the original files went

| Previous location | Current location |
|---|---|
| `kinematic tessellation optimization/function/` | `src/kinematics/` |
| `kinematic tessellation optimization/main.m` | `legacy/kinematic_main.m` |
| `kinematic tessellation optimization/results/` | `legacy/kinematic-cases/` |
| `kinematic tessellation optimization/forward_analysis/` | `extras/forward-analysis/` |
| `matlab-comsol joint optimisation/GA/*.m` | `src/comsol/` |
| `matlab-comsol joint optimisation/GGA*.m` | `legacy/comsol/` |

See `source-manifest.csv` for every source file and its SHA-256 before cleanup. `source-commit.txt` records the baseline. Generated traces, desktop metadata, editor backups and exploratory live scripts are omitted from the cleaned tree and remain in Git history. The original local working directory is unchanged.
