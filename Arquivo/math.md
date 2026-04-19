---
layout: default
title: Mathematics
nav_order: 3
---

<script>
MathJax = {
  tex: { inlineMath: [['$','$'], ['\\(','\\)']] }
};
</script>
<script src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js" async></script>

# Mathematical background

## The colony model

Fix the coordinate space $\mathbb{Z}^2$, with the $L^1$ metric
$$
d(p, q) \;=\; |p_1 - q_1| + |p_2 - q_2|.
$$

A **cell** is a finite, non-empty subset $c \subseteq \mathbb{Z}^2$. A **colony** is a finite set of cells
$$
C \;=\; \{c_1, \dots, c_n\}, \quad c_i \in \mathcal{P}_{\text{fin}}(\mathbb{Z}^2).
$$

The **expansion operator** promotes each coordinate of each cell to a singleton cell and re-aggregates:
$$
\mathrm{expand}(C) \;=\; \bigcup_{c \in C} \;\bigcup_{p \in c} \{\, \{p\} \,\}.
$$

Equivalently, in Lean:

```lean
def Colony.expand (C : Colony) : Colony :=
  ⟨C.cells.biUnion fun c => c.coords.image fun p => { coords := {p} }⟩
```

## Structural theorem

**Theorem (extensivity).** For every colony $C$,
$$
C \;\subseteq\; \mathrm{expand}(C).
$$

In Lean this is the lemma `Colony.expand_mono`, and it is the anchor of every quantitative lemma that follows.

## Stage calculus (v0.2)

Let $\mathrm{expand}^n$ denote the $n$-fold iteration. Define the **stage** at $n$ as
$$
\mathrm{stage}_n(C) \;=\; \mathrm{expand}^n(C).
$$

### Bound

Let $\nu$ be the size of the unit neighborhood (e.g. $\nu = 5$ for the von Neumann neighborhood including the origin). Then
$$
\lvert \mathrm{stage}_n(C) \rvert \;\le\; \lvert C \rvert \cdot \nu^n.
$$

In Lean:

```lean
theorem Colony.stage_bound
    (C : Colony) (n : ℕ) :
    (Colony.expand^[n] C).cells.card ≤ C.cells.card * (neighborhood_size ^ n)
```

### Coverage

Let $B_n(s) = \{p \in \mathbb{Z}^2 : d(p, s) \le n\}$ be the $L^1$-ball. If $s$ is a seed cell of $C$, then every point of $B_n(s)$ is covered by $\mathrm{stage}_n(C)$:
$$
\forall p \in B_n(s),\;\exists c \in \mathrm{stage}_n(C),\; p \in c.
$$

## The C → K → F → U operator chain

The colony model supplies a concrete testbed for the AXLE operator calculus. Each operator is an endofunctor on colonies:

| Operator | Symbol | Role |
| :--- | :---: | :--- |
| Compression | $\mathrm{C}$ | Project the coordinate footprint onto a quotient. |
| Curvature | $\mathrm{K}$ | Apply the local curvature correction at each cell. |
| Fold | $\mathrm{F}$ | Re-assemble corrected cells into a single colony. |
| Unification | $\mathrm{U}$ | Close the circuit: $\mathrm{U} = \mathrm{F} \circ \mathrm{K} \circ \mathrm{C}$. |

The **compatibility theorem** (v0.2 target) states that $\mathrm{U}$ commutes with the stage:
$$
\mathrm{U}\bigl(\mathrm{stage}_n(C)\bigr) \;=\; \mathrm{stage}_n\bigl(\mathrm{U}(C)\bigr).
$$

This is the sense in which the colony model is *orthogenetic*: expansion preserves the operator chain, so the dynamics are compatible with the structural closure.

## Why v0.1.0 matters

v0.1.0 is the first release in which the structural layer — `Cell`, `Colony`, `expand`, `expand_mono` — compiles as a single coherent Lean library. It is the minimal state from which the stage calculus and the AXLE bridge can be built with confidence.
