# Validation record

Checked on 8 October 2026 with MATLAB R2025a Update 1 on macOS. These results cover the maintained-code examples, not a rerun of every paper result.

## Automated checks

- 3/3 geometry tests passed: panel-edge invariance, the rotating-square expansion law, and node/cell conversion.
- 14/14 existing forward-analysis tests passed after updating their folder paths.
- The test suite is asserted to be nonempty, so an empty test discovery cannot count as a pass.

## Inverse-design examples

Both runs used a 4 × 4 unit grid, a circular target, the recorded default example options and feasibility tolerance `1e-4`.

| Quantity | Rigid | Non-rigid |
|---|---:|---:|
| Maximum nonlinear equality residual | 1.43e-5 | 9.13e-5 |
| Maximum positive inequality residual | 0 | 3.14e-6 |
| Maximum linear hinge residual | 2.48e-10 | 5.91e-13 |
| Solver exit flag | 2 | 2 |
| Feasible at stated tolerance | Yes | Yes |
| Approximate elapsed time on this machine | 59 s | 13 s |

Runtime is indicative, not a benchmark guarantee. MATLAB reported near-singular linear systems during parts of the optimisation; the final results satisfied the reported feasibility checks. Feasibility does not establish a global optimum or fabrication accuracy.

`demo_deployment` and the rigid example generated the new README illustrations. The figures were visually inspected. Documentation links and the cleaned file tree were checked.

## Not rerun

- Full MATLAB–COMSOL genetic optimisation and finite-element studies.
- Every legacy case, parameter sweep or fabrication export.
- Exact reproduction of all published figures and experiments.

These are explicitly separate from the verified quick start.
