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
else needs setting up. Out of the box it builds a 5×5 unit pattern, writes
`output/tessellation_initial_10X10.svg` and plots the compact and deployed
states. The optimisation call is commented out — uncomment the block under
*Kinematic optimisation* in `main.m` to run it.

### Knobs in `main.m`

| Variable | Meaning |
|---|---|
| `n_rows`, `n_cols` | number of units along y and x |
| `unit_length` | side length of one unit (= 2 quadrant squares) |
| `opening_angle` | rotation angle ξ of the uniform deployment used as the initial guess |
| `target_shape` | 1 circle · 2 ellipse · 3 vase · 4 wavy · 5 heart (see `shape.m`) |
| `shape_size` | characteristic size *r* of the target curve |
| `deployability` | 1 rigid-deployable · 2 non-rigid deployable |
| `ligament_length` | length *t* of the ligament left at each hinge (SVG export only) |
| `cut_width` | width *w* of the cut void (SVG export only) |

Requires the **Optimization Toolbox** (`fmincon`).

---

## 2. Workflow

```
       main.m
         │
         ├─►  tessellation_deployment(m,n,L,0)         compact / as-cut state
         │        └─ create_tessellation ─ create_unit ─ transform_square
         │
         ├─►  tessellation_deployment(m,n,L,ξ)         uniformly deployed state
         │                                             (initial guess)
         │
         ├─►  create_svg_tessellation(...)      ──►  output/*.svg
         │
         ├─►  tessellation_optimization(guess, s, r, p)
         │        │
         │        │  objective   1/M · (l/π · angle_diff + edge_diff)
         │        │  linear eq   Aeq  — the hinge conditions
         │        │  nonlinear   rigid.m  (p = 1)   or   nonrigid.m  (p = 2)
         │        │                 ├─ angle conditions      Eq. (5), (10)
         │        │                 ├─ edge conditions       Eq. (4), (9)
         │        │                 ├─ symmetry
         │        │                 ├─ square/compact cond.  Eq. (8)
         │        │                 ├─ non-overlap           Eq. (6)
         │        │                 └─ boundary shape        Eq. (7)  ─ shape.m
         │        └─►  fmincon (interior-point)      optimised deployed state
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

The tie list is imposed as the linear equality matrix `Aeq` in
`tessellation_optimization.m` (`index1` / `index2`). The full map is
documented in the header of `create_unit.m`.

**Flat vs. cell layout.** `fmincon` works on one flat `(m·n·16) × 2` array;
every geometric routine wants the cell array. Converting between the two is
`units_to_nodes` and `nodes_to_units`. Unit `(i,j)` occupies rows
`(i-1)·n·16 + (j-1)·16 + (1:16)`. Because `fmincon` flattens a matrix variable
column-major, the *y* coordinate of node `k` sits at column `k + m·n·16` — that
is where the `+ m*n*16` offsets in `Aeq` come from.

---

## 4. File reference

### Top level

| File | Kind | Purpose |
|---|---|---|
| `main.m` | script | driver; sets parameters and calls the pipeline |
| `Energy_rigid.m` | script | bending energy of a rigid-deployable design |
| `Energy_nonrigid.m` | script | bending energy of a non-rigid design |
| `deployment.m` | script | scratch: replay the deployment of one unit |
| `*.mat` | data | saved optimised patterns (`ellipse_rigid`, `vase_rigid`, `wavy`, `wavy2`) |
| `output/` | data | exported SVG cut patterns |

### `function/` — geometry and pipeline

| File | Purpose |
|---|---|
| `create_unit.m` | split one square into 4 quadrant squares — **defines the node numbering** |
| `create_tessellation.m` | build the m×n grid of compact units |
| `tessellation_deployment.m` | uniform deployment by opening angle ξ; also produces the compact state at ξ = 0 |
| `transform_square.m` | rotate about a pivot + translate, Eqs. (1)–(2) |
| `tessellation_optimization.m` | sets up and runs `fmincon` |
| `tessellation_compaction.m` | fold a deployed pattern back to its as-cut state |
| `create_svg_tessellation.m` | export the compact pattern as a cuttable SVG |

### `function/` — primitives

| File | Purpose |
|---|---|
| `angle_calculate.m` | unsigned corner angle at the middle of three nodes |
| `length_calculate.m` | distance between two nodes |
| `ifoverlapping.m` | signed corner orientation `⟨v₁×v₂, n̂⟩`, Eq. (6) |
| `angle_diff.m` | objective term: angle variation between neighbouring units |
| `edge_diff.m` | objective term: edge-length variation between neighbouring units |
| `nodes_to_units.m` / `units_to_nodes.m` | flat ⇄ cell layout conversion |
| `shape.m` | implicit equation of the target boundary curve |
| `plot_tessellation.m` | draw the rigid panels |
| `plot_boundary.m` | draw the target curve |

### `function/constraints/`

| File | Purpose |
|---|---|
| `rigid.m` | nonlinear constraints, rigid-deployable (parallelogram voids) |
| `nonrigid.m` | nonlinear constraints, non-rigid deployable (weaker conditions) |
| `rigid.asv` | MATLAB autosave, not used |

Both return `[c, ceq]` in the block order: angle → edge → symmetry → square
condition → non-overlap → boundary. `fmincon` enforces `c ≤ 0` and `ceq = 0`.

### `function/shapes/`

Standalone helper scripts for individual target curves (`circle`, `egg`,
`star`, `wavy`, `wedge`, `anvil`, `rainbow`, `rect`, `shear`). The pipeline
itself reads its target from `shape.m`, not from these.

### Not on the active path

`tessellation_con_notuse.m` (superseded by `rigid.m`/`nonrigid.m`),
`main.asv`, `rigid.asv`, `test.mlx`, `gif.mlx`, `case1.fig`.

---

## 5. Rigid vs. non-rigid deployability

| | non-rigid (`p = 2`) | rigid (`p = 1`) |
|---|---|---|
| intersections per unit | 1 | 2 |
| edge condition | `aᵢ − a'ᵢ = 0`, Eq. (4) | `aᵢ − bᵢ = 0`, Eq. (9) |
| angle condition | `Σθᵢ = 2π`, Eq. (5) | additionally `θ₁+θ₂ = θ₃+θ₄ = π`, Eq. (10) |
| voids | general quadrilaterals | **parallelograms** |
| deployment | panels may be geometrically frustrated | pure rotation, panels keep size and angle |

Both also impose the non-overlap (Eq. 6), boundary shape (Eq. 7) and compact
shape (Eq. 8) conditions.

---

## 6. Notes and known issues

* **`shape.m` and `plot_boundary.m` duplicate the curve parameters.** `shape.m`
  drives the optimiser, `plot_boundary.m` only draws. They are not read from a
  common source, so a curve changed in one must be changed in the other. The
  wavy curve in particular is `0.3·cos(πx) + 1.2` in `shape.m` but
  `0.45·cos(0.8πx) + 1.5` in `plot_boundary.m` and in the `s == 4` branch of the
  constraint files — the constraint branch is what actually governs the result.
* **Suspected copy-paste slip in `Energy_rigid.m` / `Energy_nonrigid.m`.** In the
  "bending hinge in adjacent units (vertical)" block the compacted triple mixes
  units `(i+1,j)` and `(i,j)` while the deployed triple is read entirely from
  `(i+1,j)`. The third deployed node is probably meant to come from `ds{i,j}`.
  Flagged in a comment in both files and deliberately left unchanged, since
  fixing it would change published energy values.
* **`Aeq` is a dense `2mn16 × 2mn16` matrix.** For a 5×5 grid that is 800×800,
  but it grows as the fourth power of the grid size. The rows past the hinge
  conditions are all zero and trivially satisfied. Switching to `sparse` would
  be the single biggest win if larger grids are ever needed.
* **`c` and `ceq` grow inside loops** in the constraint files. Harmless but the
  main remaining Code Analyzer warning; preallocating means counting the
  constraints up front.
* **The symmetry block assumes even `m` and `n`.** The loops run `1:m/2` and
  `1:n/2`, so an odd grid silently drops the middle row/column of symmetry
  conditions.
* The default `MaxFunEvals` of 8000 is often hit before convergence; the solver
  reports "stopped prematurely". Raise it in `tessellation_optimization.m` for
  production runs.
