---
layout: default
title: Roadmap
nav_order: 4
---

# Roadmap

The canonical roadmap lives at the repository root in [`ROADMAP.md`](https://github.com/<ORG>/orthogenesis/blob/main/ROADMAP.md). This page renders the same sequencing for the documentation site.

## Phase 0 — Structural layer *(done in v0.1.0)*

- ✅ `Cell.ext`
- ✅ `Colony.ext`
- ✅ `Colony.expand_mono`

## Phase 1 — Stage calculus *(v0.2)*

- 🧭 `Colony.expand_card_le`
- 🧭 `Colony.stage_bound`
- 🧭 `Cell.expand_covers_neighborhood`
- 🧭 `Colony.coord_coverage`

## Phase 2 — AXLE bridge *(v0.3)*

- 🚧 `AXLE.compatibility`
- 🚧 `AXLE.fold_unification`

## Phase 3 — Reproducibility

- CI pipeline with `lake exe cache get`.
- Pinned Mathlib revision in `lake-manifest.json`.
- `leanblueprint` rendered as `blueprint.pdf`.

## Phase 4 — Integration

- Publish `docs/` as a GitHub Pages site.
- Drop the Orthogenesis block into the Geometry Hub.
- Mint a Zenodo DOI for v0.1.0.
- SBM Bienal poster (shipped with the release archive).
