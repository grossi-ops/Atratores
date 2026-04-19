---
layout: default
title: Orthogenesis
nav_order: 1
---

<script>
MathJax = {
  tex: { inlineMath: [['$','$'], ['\\(','\\)']] }
};
</script>
<script src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js" async></script>

# Orthogenesis

*A Lean 4 formalization of the colony / cell expansion model and the C → K → F → U operator chain.*

[Source on GitHub](https://github.com/<ORG>/orthogenesis){: .btn } [Release v0.1.0](https://github.com/<ORG>/orthogenesis/releases/tag/v0.1.0){: .btn }

---

## What this is

Orthogenesis is a machine-checked development, written in Lean 4 and built on top of Mathlib, of the discrete-geometric model that underlies the **C → K → F → U** operator chain:

$$
\underbrace{\mathrm{C}}_{\text{Compression}}
\;\longrightarrow\;
\underbrace{\mathrm{K}}_{\text{Curvature}}
\;\longrightarrow\;
\underbrace{\mathrm{F}}_{\text{Fold}}
\;\longrightarrow\;
\underbrace{\mathrm{U}}_{\text{Unification}}.
$$

The project is organized around three layers:

1. **Structural** — the canonical `Cell`, the finite `Colony`, and the `expand` operator.
2. **Stage calculus** — iteration of `expand` with the quantitative lemmas `stage_bound` and `coord_coverage`.
3. **AXLE bridge** — the operator calculus connecting the colony model to C → K → F → U.

The `v0.1.0` release closes the structural layer: every definition compiles, extensionality holds, and `expand_mono` is proved.

## Key result (structural layer)

**Theorem (`Colony.expand_mono`).** For every colony $C$,

$$
C.\mathrm{cells} \;\subseteq\; (\mathrm{expand}\ C).\mathrm{cells}.
$$

That is, one step of expansion is extensive. All subsequent quantitative statements (`stage_bound`, `coord_coverage`) rest on this lemma.

## Navigation

- **[Modules](./modules.html)** — file-by-file reference of the Lean sources.
- **[Mathematical background](./math.html)** — the discrete-geometric model, rendered with MathJax.
- **[Roadmap](./roadmap.html)** — sequencing of the open lemmas for v0.2 and beyond.

## Build status

```
✔ [1731/1732] Built Orthogenesis
Build completed successfully.
```

## Citing Orthogenesis

```bibtex
@software{orthogenesis_v0_1_0,
  author  = {<AUTHOR NAME>},
  title   = {Orthogenesis: A Lean 4 formalization of the colony expansion model},
  version = {v0.1.0},
  year    = {2026},
  doi     = {<DOI once minted>},
  url     = {https://github.com/<ORG>/orthogenesis}
}
```
