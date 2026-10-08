# Workflow

## 1. Kinematic inverse design

```matlab
setup_project
cfg = kinematic_config;
cfg.target = make_target('vase', 'Size', 3.2, 'ClampX', 2.5);
rigid = main_rigid(cfg);
nonrigid = main_nonrigid(cfg);
```

Both modes use the same target and settings. The constraint files in `src/kinematics/constraints/` retain the separate deployability conditions. The result contains the initial guess, deployed geometry, compact pattern, solver settings and feasibility residuals.

Targets live in `src/shapes/`:

| Name | Boundary | Mechanical stage |
|---|---|---|
| `circle` | Radius `Size` | Not supported |
| `ellipse` | Semi-axes `Size` and `Size/2`, vertical grips at ±`ClampX` | Supported |
| `vase` | Opposed circular arcs, vertical grips | Supported |
| `wavy` | ±[0.45 cos(0.8π(x−1.25)) + 1.5] | Supported |
| `heart` | Implicit heart curve | Not supported |

`make_target` also accepts an implicit function `@(x,y)` (zero on the boundary, negative inside), or an N-by-2 closed polyline. Custom shapes are currently kinematic-only. Edit `target_profile.m` when changing an explicit upper/lower profile, and keep its implicit equation in `shape.m` consistent. The wavy profile is now shared by both stages; the old mechanical script's separate amplitude of 0.2 has been removed.

For a kinematic-only example with freely varying sheet size:

```matlab
cfg.rows = 4; cfg.cols = 4;
cfg.target = make_target('circle','Size',3.2*sqrt(2/2.2));
cfg.unit_length = [];                  % fit the initial scale
cfg.solver.FreeScale = true;
cfg.solver.Symmetry = false;
cfg.solver.FeasibilityTolerance = 1e-4;
result = main_rigid(cfg);
```

A converged geometry satisfies the reported tolerance; this does not establish a globally optimal design. Different targets may require another initial guess or grid. Infeasible output is saved for inspection but is not accepted by the mechanical stage.

## 2. MATLAB–COMSOL refinement

Start a COMSOL Multiphysics server using your installation, add its `mli` folder to the MATLAB path, and connect to the port it reports:

```matlab
addpath(fullfile(comsol_installation, 'mli'))
mphstart(server_port)
settings = mechanical_config;
mechanical = main_comsol(rigid, settings);
```

`comsol_installation` and `server_port` are your local installation path and server port. They are intentionally not stored in the repository. You can also start MATLAB through your installed COMSOL–MATLAB launcher.

The hand-off checks that the input is feasible, its compact sheet has size `cols × rows`, and its full geometry agrees with the symmetric mechanical parameterisation. It never silently replaces an asymmetric pattern by a mirrored quarter. Rigid patterns use the original 12-parameter, 2 × 4 construction; non-rigid patterns use unique quarter-grid coordinates. The GA keeps this parameterisation but does not reimpose every stage-1 boundary constraint on every candidate.

### Physical settings

- `metres_per_unit` converts kinematic coordinates to metres. The default `1` preserves the original scripts' metre scale; set this to your specimen's scale.
- `thickness` is always in metres. `cut_width`, `ligament_gap` and per-grip `displacement` use kinematic coordinate units.
- Density is kg/m³; Young's modulus is Pa; Poisson's ratio is dimensionless. The inherited example values are **not a material calibration**.
- By default, displacement is derived from `target.clampx − cols/2`. Both grips move equally and oppositely.
- The source model uses 2D solid mechanics, a hyperelastic material feature, geometric nonlinearity, pointwise horizontal grip displacement and rigid-motion suppression. Mesh size and load increments are explicit settings. It does not model out-of-plane buckling or self-contact.

Each candidate rebuilds the slit geometry with finite ligament gaps, meshes the sheet, ramps the displacement, and samples the loaded lower edge. The fitness is the area between its polynomial fit and the target profile after vertical alignment at the left grip, as in the source objective. Grip points are excluded, so the polynomial extrapolates to the grips; inspect the saved sample extent and fitted curve when assessing a result. Error is measured in squared design-coordinate units. Numerical integration replaces the old Symbolic Toolbox calculation.

### GA settings and outputs

`settings.ga` controls population, generations, perturbation radius, crossover, mutation, shape-error tolerance and stall limit. A serial real-valued GA uses tournament selection and one retained elite. This replaces the separate legacy selection/crossover/mutation drivers; it is not an identical random trajectory to the old roulette-selection code.

The first candidate is the supplied kinematic design. If this baseline cannot solve, the run stops with the COMSOL error so that material, mesh, geometry and loading can be checked. Later invalid candidates receive infinite fitness and their errors are logged; they are not silently substituted by the current best design.

Each run gets a new directory under `results/` with `inputs.mat`, `checkpoint.mat` and `mechanical_result.mat`. Results contain the best chromosome, reconstructed nodes, sampled and fitted boundary, fitness history, failure log, settings and stopping reason. A stopped run is not necessarily within the target tolerance.

## Checks

```matlab
setup_project
assertSuccess(runtests('tests'))
```

The fast tests do not require COMSOL. Actual kinematic and COMSOL runs are recorded separately in `validation.md`. Keep generated plots, datasets and model files in `results/`, not in source folders.
