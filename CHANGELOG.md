# Changes

## Focused two-stage project — October 2026

- Three entry points: rigid kinematics, non-rigid kinematics, and COMSOL/GA refinement.
- Removed unused forward studies, duplicate drivers, saved legacy workspaces and the large solver reference.
- Consolidated named target shapes and fixed built-in target structs; shared the wavy curve between stages.
- Preserved the kinematic constraint equations and the original mechanical cut constructions.
- Replaced workspace-dependent GA scripts with a bounded serial driver, explicit settings, checkpoints and visible failure handling.
- Removed material/loading overrides, corrected the rigid mirror function name and non-rigid grid indexing, and used reference coordinates plus displacement for loaded boundary extraction.
- Removed unused bending-energy extraction and per-candidate plotting. COMSOL models are rebuilt from MATLAB; `.mph` files are ignored.

Based on maintained research source commit `a567cb4` and the first cleanup commit `632894c`. Original scripts remain in Git history and in the untouched research folder. This is a maintained version, not a frozen paper reproduction.
