# Kinematic tessellation optimisation

MATLAB implementation of **stage 1** of the inverse-design framework in

> X. Ying, D. Fernando, M. A. Dias, *Inverse design of programmable shape-morphing
> kirigami structures*, International Journal of Mechanical Sciences **286** (2025) 109840.
> <https://doi.org/10.1016/j.ijmecsci.2024.109840>

Given a target boundary shape, the code finds a quad-kirigami cut pattern that
deploys onto that shape while remaining geometrically valid (Section 2.1,
Fig. 1 and Fig. 2 of the paper). The result is the initial cut pattern for
stage 2, the GA + FEA loop that adds elasticity.

---

## 1. Quick start

```matlab
cd '.../kinematic tessellation optimization'
main
```

`main.m` adds `function/` (and its sub-folders) to the path itself, so nothing
else needs setting up. Requires the **Optimization Toolbox** (`fmincon`).

**Before a long run, check the problem is solvable at all:**

```matlab
design_freedom(4, 4, 1, 3.2*sqrt(2.0/2.2), 1)
```

If that reports zero remaining freedom, no amount of solver effort will help —
see [section 6](#6-making-the-optimiser-converge).

### Knobs in `main.m`

| Variable | Meaning |
|---|---|
| `n_rows`, `n_cols` | number of units along y and x. **Keep both even** — the symmetry conditions loop over `1:m/2`, so an odd grid silently drops the middle row and column |
| `target_shape` | 1 circle · 2 ellipse · 3 vase · 4 wavy · 5 heart (see `shape.m`) |
| `shape_size` | characteristic size *r* of the target curve |
| `deployability` | 1 rigid-deployable · 2 non-rigid deployable |
| `FIT_UNIT_LENGTH` | `true` fits the unit length as well as the opening angle, and lets the optimiser keep tuning the overall scale (`FreeScale`) |
| `unit_length` | the fixed size when `FIT_UNIT_LENGTH` is `false`, otherwise the starting size |
| `ligament_length` | length *t* of the ligament left at each hinge (SVG export only) |
| `cut_width` | width *w* of the cut void (SVG export only) |

---

## 2. Workflow

```
       main.m
         │
         ├─►  design_freedom(...)                  IS THIS SOLVABLE?  (seconds)
         │
         ├─►  fit_initial_guess(m,n,s,r,mode)      pick opening angle + size
         │        └─ boundary_residual ─ shape     closest to the target
         │
         ├─►  tessellation_deployment(m,n,L,0)     compact / as-cut state
         │        └─ create_tessellation ─ create_unit ─ transform_square
         │
         ├─►  create_svg_tessellation(...)      ──►  output/*.svg
         │
         ├─►  tessellation_optimization(guess, s, r, p, opts)
         │        │
         │        │  objective   1/M · (l/π · angle_diff + edge_diff)
         │        │  linear eq   Aeq  — the hinge conditions (sparse)
         │        │  nonlinear   rigid.m  (p = 1)   or   nonrigid.m  (p = 2)
         │        │                 ├─ angle conditions      Eq. (5), (10)
         │        │                 ├─ edge conditions       Eq. (4), (9)
         │        │                 ├─ symmetry              OPTIONAL, off
         │        │                 ├─ square/compact cond.  Eq. (8)
         │        │                 ├─ non-overlap           Eq. (6)
         │        │                 ├─ boundary shape        Eq. (7)
         │        │                 │      └─ boundary_residual
         │        │                 └─ min edge length       optional
         │        └─►  fmincon (interior-point) + warm restarts
         │                  └─► info.converged   ← ALWAYS CHECK THIS
         │
         ├─►  tessellation_compaction(optimised)       fold back to as-cut state
         │
         └─►  plot_tessellation + plot_boundary        figures
```

Downstream, once you have both states in the workspace:

* `Energy_rigid.m` / `Energy_nonrigid.m` — rotating-spring bending energy, Eqs. (11), (14)
* `deployment.m` — scratch script that replays the motion of one unit

---

## 3. The 16-node unit — read this first

Everything in the code is index arithmetic on one data structure: the
tessellation is an `m × n` **cell array**, and each cell holds a **16 × 2**
list of node coordinates. Rows 1:4, 5:8, 9:12 and 13:16 are the four rigid
quadrant squares Q1…Q4 of one unit, each listed clockwise.

```
     14 --------15/2-------- 3
      |    Q4    |    Q1     |          Q1 = rows 1:4   top-right
      |          |           |          Q2 = rows 5:8   bottom-right
     13        16/1          4          Q3 = rows 9:12  bottom-left
     /10        11/6         7\         Q4 = rows 13:16 top-left
      |    Q3    |    Q2     |
      |          |           |
      9 --------12/5-------- 8
```

Nodes shared by two squares are **stored twice**, once per square. Which
duplicates are tied together is exactly what makes this a kirigami:

| Pair | Meaning |
|---|---|
| `4↔7`, `10↔13` | mid-edge hinges inside a unit |
| `1↔16`, `6↔11` | the two centre hinges inside a unit |
| `8↔9`, `3↔14` | hinges to the unit on the right, `(i,j)→(i,j+1)` |
| `15↔12`, `2↔5` | hinges to the unit above, `(i,j)→(i+1,j)` |
| `2↔15`, `5↔12` *within* a unit, and the two centre pairs relative to each other | **not tied** — these open into the rotating-squares void |

Boundary nodes, per unit: left `14, 9` · right `3, 8` · bottom `12, 5` · top `15, 2`.

The tie list is imposed as the linear equality matrix `Aeq` in
`tessellation_optimization.m` (`index1` / `index2`). The full map is in the
header of `create_unit.m`.

**Flat vs. cell layout.** `fmincon` works on one flat `(m·n·16) × 2` array;
every geometric routine wants the cell array. `units_to_nodes` and
`nodes_to_units` convert between them. Unit `(i,j)` occupies rows
`(i-1)·n·16 + (j-1)·16 + (1:16)`. Because `fmincon` flattens a matrix variable
column-major, the *y* coordinate of node `k` sits at column `k + m·n·16` — that
is where the `+ m*n*16` offsets in `Aeq` come from.

---

## 4. Constraint blocks

`rigid.m` and `nonrigid.m` share the same structure and both return
`[c, ceq]`; `fmincon` enforces `c ≤ 0` and `ceq = 0`. Both take an optional
6th argument, a struct of options.

| Block | Condition | Paper | Kind |
|---|---|---|---|
| 1 · angle | Angles around every interior corner sum to 2π. In `rigid.m`, pairs across a cut additionally sum to π (colinear when closed). | 5, 10 | `ceq` |
| 2 · edge | Corresponding edges either side of a cut have equal length; with block 1 this makes every void a parallelogram in the rigid case. | 4, 9 | `ceq` |
| 3 · symmetry | Left/right and top/bottom mirror symmetry. **Optional, off by default** — see section 6. | — | `ceq` |
| 4 · square | The compacted state closes into an n×m rectangle. With `FreeScale`, only the aspect ratio is pinned. | 8 | `ceq` |
| 5 · non-overlap | `⟨v₁×v₂, n̂⟩ ≥ 0` at every corner, plus all 16 corners of each quadrant square. | 6 | `c` |
| 6 · boundary | Boundary nodes on the target curve, interior nodes inside it. Delegated to `boundary_residual.m`. | 7 | `c`, `ceq` |
| 7 · min edge | Optional floor on every panel edge — see section 6. | — | `c` |

### Options

| Option | Default | Effect |
|---|---|---|
| `Symmetry` | `false` | Add the mirror conditions. Over-constrains the problem; see section 6. |
| `FreeScale` | `false` | Replace the two absolute-size equalities with one aspect-ratio equality, freeing the overall scale. |
| `MinEdgeLength` | `0` (constraints) / `0.15 × unit` (via `tessellation_optimization`) | Lower bound on every quadrant-square edge. |

---

## 5. Rigid vs. non-rigid

| | non-rigid (`p = 2`) | rigid (`p = 1`) |
|---|---|---|
| intersections per unit | 1 | 2 |
| edge condition | `aᵢ − a'ᵢ = 0`, Eq. (4) | `aᵢ − bᵢ = 0`, Eq. (9) |
| angle condition | `Σθᵢ = 2π`, Eq. (5) | additionally `θ₁+θ₂ = θ₃+θ₄ = π`, Eq. (10) |
| voids | general quadrilaterals | **parallelograms** |
| deployment | panels may be geometrically frustrated | pure rotation, panels keep size and angle |

Both also impose non-overlap (6), boundary shape (7) and compact shape (8).

---

## 6. Making the optimiser converge

The optimiser used to return badly distorted patterns whose compacted state was
not a rectangle. That had **three separate causes**, in order of importance.

### 6.1 The problem was over-constrained (the real one)

The design variables are the `2·m·n·16` node coordinates. Every *independent*
equality removes one dimension. Counting the rank of the equality Jacobian for
a 4×4 circle in **rigid** mode, with the symmetry conditions on:

| block | rows | rank added | freedom left |
|---|---|---|---|
| hinges (linear) | 224 | 224 | 288 |
| angle conditions | 93 | 93 | 195 |
| edge conditions | 112 | 112 | 83 |
| symmetry left-right | 128 | 56 | 27 |
| symmetry top-bottom | 64 | 19 | 8 |
| square conditions | 36 | 8 | **0** |
| **boundary target curve** | **32** | **0** | **0** |

By the time the target-shape equations are reached there is **nothing left to
move**, so they add no rank and cannot be satisfied. 689 equations for 512
unknowns. The solver was not slow — the problem had no solution, and what came
back was a least-squares compromise between incompatible demands. That is
exactly what collapsed panels and a non-rectangular compacted state look like.

The symmetry conditions are the cause. They are not a physical requirement —
the paper introduces them only to cut the parameter count — and more than half
their rows are redundant with each other anyway. **They are now off by
default**, which restores:

| mode | symmetry on | symmetry off | off + `FreeScale` |
|---|---|---|---|
| rigid 4×4 | 0 | **19** | **20** |
| non-rigid 4×4 | 0 | **59** | **60** |

Run `design_freedom(m, n, s, r, p)` on your own case before a long solve. If
you genuinely need a symmetric design, the sound way is to optimise a quarter
of the sheet and mirror it — reducing the *variables* rather than adding
equations, which is what the paper does for the FEA stage.

### 6.2 The evaluation budget was far too small

`MaxFunEvals` was hard-coded to 8000. `fmincon` builds gradients by finite
differences, so one iteration costs `nvars + 1` evaluations — for a 4×4 grid
(512 variables) that is **about fifteen iterations** before the solver gives up.
The default is now `400 × nvars`, and `Restarts` (default 3) warm-restarts the
solver if it still hits a limit. On the non-rigid 4×4 circle, budget alone took
`max|ceq|` from 9.0e-03 to 2.6e-06 and the panel area ratio from 18 to 4.8.

`Aeq` was also allocated dense at `2mn16 × 2mn16` with only a fraction of the
rows used — 224 real rows and 1824 all-zero ones for a 4×4 grid. It is now
built sparse and trimmed, which removes that rank deficiency and the
`O(grid⁴)` memory cost.

### 6.3 `angle_calculate` could return complex numbers

The law-of-cosines form, `acos((a²+b²−c²)/(2ab))`, can have its argument pushed
just outside `[-1,1]` by rounding, and MATLAB's `acos` then returns a **complex
number** that propagates straight into the constraint vector. Probing
near-degenerate corners — exactly the flat and folded configurations an
optimiser passes through — **6019 of 400000 triples came back complex**, and
coincident nodes returned `NaN`. Its derivative is also unbounded there, so
finite-difference gradients at those corners are meaningless.

It now uses `atan2(|u×v|, u·v)`, which agrees to ~1e-11 on well-conditioned
corners, is well conditioned over the whole range, and returns 0 rather than
`NaN` for coincident nodes.

### 6.4 Other fixes

* **Panel collapse.** The non-overlap conditions only fix the *sense* of each
  corner, so a panel can shrink towards zero area — the rigid 4×4 circle had a
  400:1 ratio between its largest and smallest panels. `MinEdgeLength` puts a
  floor under every edge. It defaults to `0.15 × unit size`; raise it if panels
  still collapse, set it to 0 if it distorts the result.
* **Objective coverage.** `angle_diff` and `edge_diff` looped `i = 1:m-1,
  j = 1:n-1` and then patched in a single extra term for unit `(m,n)`, which
  left most of the last row and column out of the objective — precisely the
  units carrying the boundary condition, free to distort at no cost. Both now
  cover every adjacent pair exactly once. **This changes the objective value**,
  so numbers will not match the published run exactly.

### 6.5 Tuning the size to fit the boundary

Section 2.1.2 of the paper says the opening angle is chosen by "minimising the
difference between the initial configuration and the target shape". The code
never did this — the angle was hard-coded. `fit_initial_guess` now does it, and
extends the same idea to the unit length:

```matlab
[guess, info] = fit_initial_guess(m, n, s, r, mode);       % fit angle and size
[guess, info] = fit_initial_guess(m, n, s, r, mode, 1.0);  % fit angle only
```

Initial boundary mismatch (RMS) on a 4×4 grid:

| shape | hand-picked π/6 | fitted angle | fitted angle + length |
|---|---|---|---|
| circle | 0.346 | 0.328 | 0.323 |
| ellipse | 1.757 | 0.896 | 0.610 |
| vase (rigid) | 1.637 | 1.339 | 0.786 |
| wavy | 0.942 | 0.499 | 0.477 |
| heart | 575.8 | 102.8 | **2.005** |

Fitting the length as well as the angle is what does the work on the harder
shapes. `fit_initial_guess` warns when the best angle lands on the search
bound — for the ellipse, vase and wavy targets it does, which means a uniform
grid is simply a poor match for those outlines and the optimiser has to do the
rest.

To let the optimiser keep adjusting the scale afterwards, pair it with
`FreeScale`, which replaces the two absolute-size equalities with a single
aspect-ratio one. `main.m` wires both to the `FIT_UNIT_LENGTH` switch.

The search deliberately excludes the ends of the mechanism's travel
(default ξ ∈ [5.7°, 84.3°]). Minimising boundary mismatch alone often drives ξ
to 0 or π/2, but there every void has zero area, so the non-overlap
inequalities sit exactly at equality — an interior-point solver started on its
own constraint boundary has nowhere to step.

---

## 7. Any boundary shape

The five hard-coded shapes are no longer a limit. `make_target` accepts three
forms, and everything downstream — `rigid.m`, `nonrigid.m`,
`fit_initial_guess`, `design_freedom`, `plot_boundary` — takes the result
wherever it used to take a shape code.

```matlab
% 1. the built-ins, unchanged
t = make_target(1, 'Size', 3.05);

% 2. any implicit curve: 0 on the boundary, negative inside
t = make_target(@(x,y) (abs(x)/3.0).^4 + (abs(y)/2.4).^4 - 1);   % superellipse
t = make_target(@(x,y) ((x.^2+y.^2).^2 - 6*(x.^2-y.^2))/10 - 1); % peanut

% 3. any closed curve given as points — traced, sampled or exported from CAD
th = linspace(0, 2*pi, 400).';  th(end) = [];
r  = 2.7*(1 + 0.18*cos(5*th));
t  = make_target([r.*cos(th), r.*sin(th)]);                      % flower

% add tensile grips: left/right held straight instead of following the curve
t = make_target(@(x,y) x.^2/9 + y.^2/4 - 1, 'ClampX', 2.5);
```

Then simply:

```matlab
target = make_target(@(x,y) (abs(x)/3).^4 + (abs(y)/2.4).^4 - 1);
design_freedom(4, 4, target, 1, 2)                  % solvable?
guess = fit_initial_guess(4, 4, target, 1, 'nonrigid');
T = tessellation_optimization(guess, target, 1, 2, struct('FreeScale',true));
```

**How it works.** Eq. (7) of the paper is already general — it asks that each
boundary node sit on the target, measured by the distance to its projection.
`target_value.m` provides exactly that one number per point: zero on the
boundary, negative inside, positive outside. It serves as the equality for
boundary nodes and the inequality for every other node. For an implicit target
it is the curve's own function; for a point list it is the signed distance to
the polyline, which *is* the projection distance.

**Fidelity.** Handing the general path the circle as an implicit function
reproduces the built-in code exactly (`max diff 0.000e+00` on both `c` and
`ceq`). The five original codes keep their hard-coded treatment, so published
results still reproduce bit-for-bit.

**Caveat on point lists.** The signed distance to a polyline is continuous but
its derivative jumps where the nearest segment changes. Finite differences
tolerate it, but sample the curve densely — 300–400 points is comfortable —
and prefer the implicit form when you have a formula. A 400-point polygon
inscribing a circle biases the residual by ~1e-4 relative to the exact circle,
which is the polygon approximation, not an error.

---

## 8. Run time — the setting that matters

**A 4×4 circle solves in ~17 s (non-rigid) or ~55 s (rigid).** If yours takes
minutes, it is almost certainly running past the point of usefulness.

`fmincon` terminates on its **step** tolerance, not on feasibility. On a 4×4
circle the geometry is solved early and everything after that is polishing the
objective — making panels more uniform:

| reach \|ceq\| < | iteration | time |
|---|---|---|
| 1e-04 | 25 | 13 s |
| 1e-06 | 33 | 16 s |
| 1e-07 | 265 | 140 s |
| stops | 360 | 195 s |

`EarlyStop` (default **true**) leaves as soon as the design is *feasible* and
the objective has flattened (<1% over 15 iterations):

| mode | EarlyStop | time | iters | max \|ceq\| | objective | panel ratio |
|---|---|---|---|---|---|---|
| non-rigid | **on** | **17 s** | 33 | 7.2e-07 | 1.0498 | 5.0 |
| non-rigid | off | 197 s | 360 | 5.5e-10 | 0.8163 | 4.6 |
| rigid | on | 162 s | 298 | 2.4e-06 | 0.9768 | 4.9 |
| rigid | off | 166 s | 298 | 2.4e-06 | 0.9768 | 4.9 |

So non-rigid gets **11× faster** for panels ~9% less uniform.

### Rigid needs a reachable target

Rigid never gets below 1e-6 — its best is 2.4e-06 — so the `EarlyStop`
feasibility gate never opens at the default and it runs to completion. That is
correct behaviour (it will not stop at a point it has not verified), but it
means rigid needs `FeasibilityTolerance` set to something it can hit:

| FeasibilityTolerance | time | max \|ceq\| | objective |
|---|---|---|---|
| 1e-03 | 10 s | 6.6e-04 | 1.2155 (visibly irregular) |
| **1e-04** | **55 s** | 1.4e-05 | **0.9780** |
| 1e-05 | 66 s | 9.8e-06 | 0.9773 |
| 1e-06 | 165 s | 2.4e-06 | 0.9768 |

```matlab
opts.FeasibilityTolerance = 1e-4;   % rigid, under a minute
```

At 1e-4 the objective is within 0.1% of the 165-second answer and the geometry
is visually identical. In model units the sheet is ~4 across, so 1e-4 is about
2.5 µm on a 100 mm sheet — roughly 40× finer than a laser kerf.

### A fixed iteration cap is not a substitute

Interior-point feasibility is not monotone. Capping at 100 iterations gave
`|ceq| = 2.6e-05` — *worse* than capping at 60 (1.2e-07). `EarlyStop` only ever
stops at a point it has checked.

---

## 9. Performance — what actually helps

A converged 4×4 solve takes **10–20 minutes**. The old 9-second runs were not
fast, they stopped after about fifteen iterations and returned a point that
satisfied nothing.

Measured cost of one evaluation, and what a solve therefore costs:

| grid | variables | constraint eval | objective eval | one gradient | ~300 iterations |
|---|---|---|---|---|---|
| 2×2 | 128 | 3.4 ms | 2.7 ms | 0.8 s | ~4 min |
| 4×4 | 512 | 6.0 ms | 1.4 ms | 3.8 s | ~19 min |
| 6×6 | 1152 | 2.6 ms | 0.9 ms | 4.0 s | ~20 min |

There are no analytic gradients, so `fmincon` finite-differences everything:
**one gradient costs `nvars + 1` full evaluations of the objective and the
whole constraint set.** That single factor is the entire performance story.

Four things that sound like they should help, measured:

| idea | measured effect |
|---|---|
| `UseParallel` over 15 cores | **0.6× — i.e. 1.6× slower.** One evaluation is only ~6 ms, too little work to cover `parfor`'s dispatch overhead. Off by default. |
| **sparse-colouring gradients** | **evaluations fall 20–40×, wall clock does not improve.** See below. Off by default (`SparseGradients`). |
| vectorising the leaf functions | ~1.7× on the leaf call itself, ≈1.3× overall. Not worth rewriting the constraint files for. |
| preallocating `c` / `ceq` | ~10% overall. |

### The sparse-gradient attempt, and why it does not pay off

The constraint Jacobian is 1.1% dense at 4×4 — each condition touches only ~6
of the 512 variables. Columns that never share a constraint row can be
perturbed *together* in one evaluation, so a graph colouring needs 73
evaluations where fmincon uses 513 (109 vs 1153 at 6×6). `sparse_fd.m`
implements this, and it is verified correct: zero pattern entries missed, and
with a matched step it reproduces the dense Jacobian to `0.000e+00`.

The saving in evaluations is real and large. The saving in **time** is not:

| variant | evaluations | iterations | wall clock | converged | objective |
|---|---|---|---|---|---|
| dense (fmincon default) | 186,755 | 360 | 191 s | yes | 0.8163 |
| sparse, forward differences | 1,765 | 318 | **94 s** | **no** (2.2e-07) | 0.8710 |
| sparse, central differences | 9,335 | 3,104 | 1024 s | yes | **0.6493** |

Supplying gradients changes fmincon's search path, and the iteration count
swings far more than the per-iteration saving is worth. Forward differences
are twice as fast but stall above the feasibility tolerance; central
differences converge but take five times as long.

One genuine side effect: **the sparse-gradient runs consistently found better
optima** (0.649 and 0.682 against 0.816). If solution quality matters more
than run time, `SparseGradients = true` with central differences is worth a
try. Untested hypothesis for the iteration blow-up: `ScaleProblem` may
interact badly with user-supplied gradients.

Beware the MATLAB profiler here: it inflates per-call cost for very small
functions, and suggested `angle_calculate` was 33% of an evaluation when
unprofiled timing puts it nearer 7%. Trust `tic`/`toc`, not `profile`.

What might still give an order of magnitude:

* **Analytic constraint gradients** (`SpecifyConstraintGradient`). This removes
  the `nvars + 1` factor outright. It is real work — roughly 900 constraints to
  differentiate — but it is the only change that turns 20 minutes into
  seconds.
* **Fewer variables.** The hinge equalities mean a 4×4 grid has only 144
  distinct points (288 coordinates) hiding behind 512 variables. Optimising the
  distinct points directly would cut gradient cost ~1.8× *and* remove all 224
  linear equalities, improving conditioning at the same time.

Cheap levers meanwhile: lower `MaxFunctionEvaluations` (the default
`400 × nvars` is deliberately generous, and `Restarts` covers hitting the cap),
loosen `FeasibilityTolerance` from 1e-7, or start from a better
`fit_initial_guess` so fewer iterations are needed.

---

## 10. File reference

### Top level

| File | Kind | Purpose |
|---|---|---|
| `main.m` | script | driver; sets parameters and calls the pipeline |
| `Energy_rigid.m` | script | bending energy of a rigid-deployable design |
| `Energy_nonrigid.m` | script | bending energy of a non-rigid design |
| `deployment.m` | script | scratch: replay the deployment of one unit |
| `*.mat` | data | saved patterns — `ellipse_rigid`, `vase_rigid`, `wavy`, `wavy2` |
| `output/` | data | exported SVG cut patterns |

### `function/` — pipeline

| File | Purpose |
|---|---|
| `create_unit.m` | splits one square into 4 quadrant squares — **defines the node numbering** |
| `create_tessellation.m` | builds the m×n grid of compact units |
| `tessellation_deployment.m` | uniform deployment by ξ; also the compact state at ξ = 0 |
| `fit_initial_guess.m` | fits the opening angle and unit length to the target |
| `tessellation_optimization.m` | sets up and runs `fmincon`; returns `info` |
| `tessellation_compaction.m` | folds a deployed pattern back to its as-cut state |
| `create_svg_tessellation.m` | exports the compact pattern as a cuttable SVG |
| `design_freedom.m` | **diagnostic**: is the problem over-constrained? |

### `function/` — primitives

| File | Purpose |
|---|---|
| `angle_calculate.m` | unsigned corner angle, `atan2` form |
| `length_calculate.m` | distance between two nodes |
| `ifoverlapping.m` | signed corner orientation `⟨v₁×v₂, n̂⟩` — Eq. (6) |
| `angle_diff.m` / `edge_diff.m` | objective terms: variation between neighbouring units |
| `nodes_to_units.m` / `units_to_nodes.m` | flat ⇄ cell layout conversion |
| `shape.m` | implicit equation of the target boundary curve |
| `boundary_residual.m` | the boundary conditions, shared by the constraints and the fitter |
| `plot_tessellation.m` / `plot_boundary.m` | drawing |

### `function/constraints/` and `function/shapes/`

`rigid.m` and `nonrigid.m` hold the nonlinear constraints; `rigid.asv` is a
MATLAB autosave. `shapes/` holds standalone helper scripts for individual
curves — the pipeline reads its target from `shape.m`, not from these.

Not on the active path: `tessellation_con_notuse.m`, `main.asv`, `test.mlx`,
`gif.mlx`, `case1.fig`.

---

## 11. Remaining issues

* **`shape.m` and `plot_boundary.m` duplicate the curve parameters.** The
  boundary *constraints* now come from `boundary_residual.m`, but
  `plot_boundary.m` still has its own hard-coded copy for drawing. The wavy
  curve is `0.3·cos(πx) + 1.2` in `shape.m` but `0.45·cos(0.8πx) + 1.5` in
  `boundary_residual.m` and `plot_boundary.m` — the latter is what governs the
  result, and `shape.m`'s `s == 4` branch is effectively dead.
* **Suspected copy-paste slip in the energy scripts.** In the "bending hinge in
  adjacent units (vertical)" block the compacted triple mixes units `(i+1,j)`
  and `(i,j)` while the deployed triple is read entirely from `(i+1,j)`. Left
  unchanged because fixing it would change published energy values.
* **`c` and `ceq` grow inside loops** in the constraint files. Harmless, but it
  is the hot path — the constraint function is called hundreds of thousands of
  times, so preallocating would be the next worthwhile speedup.
* **The symmetry block assumes even `m` and `n`** when enabled.
* **`design_freedom` costs `O((mn)²)`** because it builds the Jacobian by
  finite differences. Fine to about 6×6.
