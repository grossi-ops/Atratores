---
layout: default
title: Modules
nav_order: 2
---

# Module reference

The Lean sources are organized around a thin dependency chain so that each module can be verified in isolation.

```
Coordinates.lean
      │
      ▼
   Cell.lean
      │
      ▼
  Colony.lean
      │
      ▼
   Stage.lean        AXLE/Operators.lean  (WIP)
```

## `Orthogenesis.Coordinates`

Defines the coordinate primitive used by every higher layer.

- `abbrev Vec2 := Int × Int`
- `instance : Repr Vec2` — derived through the pair instance; replaces the previous unsafe `Real`-valued `Repr`.
- A handful of coordinate utilities (`neighbors`, `ballL1`) kept `pure` and total.

## `Orthogenesis.Cell`

The canonical definition of a single cell. Previously the source of the v0.1.0 build break when an accidental shadow copy was introduced inside `Colony.lean`.

```lean
structure Cell where
  coords : Finset Vec2
  deriving Repr, DecidableEq
```

Exposes:

- `Cell.ext` — extensionality through `coords`.
- `Cell.card` — convenience shortcut for `coords.card`.

## `Orthogenesis.Colony`

The load-bearing module. Uses `Cell` exclusively through its canonical definition — no local redefinition.

```lean
structure Colony where
  cells : Finset Cell
  deriving Repr, DecidableEq
```

Operator:

```lean
def Colony.expand (C : Colony) : Colony :=
  ⟨C.cells.biUnion fun c => (c.coords.image fun p => { coords := {p} } : Finset Cell)⟩
```

Lemmas:

- `Colony.ext`
- `Colony.expand_mono` — `C.cells ⊆ (expand C).cells`.

## `Orthogenesis.Stage`

The indexed iteration of `expand`. Contains the open statements that motivate v0.2:

- `Colony.stage_bound` (drafted).
- `Colony.coord_coverage` (drafted).

Both are sequenced in the [roadmap](./roadmap.html).

## `Orthogenesis.AXLE.Operators`

WIP scaffolding for the AXLE operator calculus. The current file declares the types of `C`, `K`, `F`, `U` as endofunctors on `Colony`, with their compatibility lemmas left as `sorry` until the stage calculus closes.

## Tests

`Tests/OrthogenesisTests.lean` is the regression smoke test that surfaces the class of failures responsible for v0.1.0:

- duplicate `Cell` ambiguity (caught by the fully-qualified `#check`);
- `Finset.insert` namespace resolution (caught by a `decide` on `card`);
- extensionality drift (caught by an explicit `example`).
