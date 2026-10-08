# Validation — 8 October 2026

Checked with MATLAB R2025a Update 1 and COMSOL 6.2, build 290, on macOS.

## Kinematic examples

Both default 2 × 4 symmetric vase examples were run with seed 1 and feasibility tolerance `1e-6`.

| Check | Rigid | Non-rigid |
|---|---:|---:|
| Maximum nonlinear equality residual | 1.95 × 10⁻⁷ | 1.07 × 10⁻⁸ |
| Maximum positive inequality residual | 0 | 0 |
| Maximum linear hinge residual | 2.81 × 10⁻¹⁰ | 2.70 × 10⁻¹² |
| Mechanical reconstruction error | 4.72 × 10⁻⁷ | 1.62 × 10⁻⁸ |
| Feasible at the stated tolerance | Yes | Yes |

The deployed and compact patterns were visually inspected. MATLAB reported near-singular intermediate systems for some runs; acceptance was based on the final residuals. Other targets and initial guesses can fail. For example, the tested ellipse starting point failed the rigid tolerance, and a non-rigid ellipse encountered a nonlinear COMSOL convergence failure. Failures are reported, not treated as successful results.

## COMSOL and GA integration

Both default vase patterns completed a short, real COMSOL-backed GA run: **population 2, one generation, perturbation radius 0.002, seed 1, stopping tolerance 0**. Other material, geometry, mesh and loading settings used `mechanical_config` defaults. Each solve ramped both grips to ±0.5 over 50 increments.

| Result | Rigid | Non-rigid |
|---|---:|---:|
| Final shape-error objective | 0.636708444 | 1.383015540 |
| Rejected candidates | 0 | 0 |
| Stopping reason | Generation limit | Generation limit |

These checks exercised model construction, meshing, full-load nonlinear solutions, loaded boundary extraction, fitness evaluation, GA iteration and saved outputs. The final geometry of each mode was separately checked for a single connected domain. Models were released from the server after evaluation; no `.mph` files were saved or included in the repository.

**This verifies the workflow, not convergence to the target shape.** The short runs did not reach the production shape-error tolerance. A full population/generation study, mesh-convergence study and specimen-specific material calibration remain separate research tasks.

## Code checks

- All 12 fast tests passed: geometry invariants, named targets, shared wavy profile, mechanical encoding, rectangular-grid indexing, invalid-input handling, shape-error translation, GA bounds/elitism and zero-fitness handling.
- MATLAB dependency analysis found no source `.m` files outside the three main entry points' dependency closure.
- MATLAB syntax checks and Git whitespace checks passed.

The numerical model uses the source COMSOL API workflow; geometry methods were checked against the [COMSOL 6.2 API reference](https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/api/com/comsol/model/GeomSequence.html).
