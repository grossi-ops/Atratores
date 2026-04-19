# 3M — Helical Attractors on Contact 3-Manifolds WEBSITE grossi-ops,github.io/Atratores/index.html ENGLISH
VEJA PORTUGUES _ ARQUIVOS PARA DOWNLOAD EM PDF E HTML VIDE REPOSITORIO

**Chapter 10 of *Principia Orthogona*, Volume IV**
**Companion site for the XII Bienal da SBM 2026 mini-curso**

Pablo Nogueira Grossi · G6 LLC · Newark, NJ · 2026
ORCID: [0009-0000-6496-2186](https://orcid.org/0009-0000-6496-2186)
Live site: <https://totogt.github.io/3M/>
Chapter DOI: [10.5281/zenodo.19379385](https://doi.org/10.5281/zenodo.19379385)

This repository hosts Chapter 10 of Volume IV of *Principia Orthogona* — *Helical Attractors on Contact 3-Manifolds* — together with the three-session **XII Bienal da SBM 2026** mini-curso delivered at UFRN (Natal-RN) and the submission artefacts Pablo Nogueira Grossi files to the Bienal in four categories. Every page is a standalone teaching artefact: the three sessions can be read in sequence as a 180-minute course, the 3D simulator runs in the browser with no install, and the DOP853 numerical lab reproduces the tables from the chapter.

---

## Citation

```bibtex
@incollection{grossi2026gtct,
  author    = {Grossi, Pablo N.},
  title     = {Helical Attractors on Contact 3-Manifolds},
  booktitle = {Principia Orthogona},
  volume    = {IV},
  year      = {2026},
  publisher = {G6 LLC}
}
```

Open-access deposit of this chapter: [*The DM3 Operator — Explicit Toy Model and Global Dynamic Analysis*](https://doi.org/10.5281/zenodo.19379385) (Zenodo, 2026).

> The DOI `10.5281/zenodo.19117400` that appears in earlier drafts of this README and in the docstring of `labs/dm3_numeric.py` is **Volume I**'s Zenodo deposit, not this chapter's. When editing `dm3_numeric.py`, update line 7 to `DOI: 10.5281/zenodo.19379385`.

---

## Live pages

| Path | What it is |
| --- | --- |
| [`index.html`](./index.html) | Publisher landing — 3D hero, "Read, Simulate, Prove, Teach" tile grid, and the Vol IV placement map |
| [`chE-gtct.html`](./chE-gtct.html) | Chapter E — *GTCT for Everyone* (nine axioms, twelve operators, four theorems, five properties of time) — warm-up primer inherited from Book 3 (Mini-Beast) |
| [`sessions/session1-contact-geometry.html`](./sessions/session1-contact-geometry.html) | Mini-curso **S1** · 60 min · Contact geometry and the dm³ system |
| [`sessions/session2-theorem-basin.html`](./sessions/session2-theorem-basin.html) | Mini-curso **S2** · 60 min · Theorem 2.1 and the asymmetric basin |
| [`sessions/session3-lean-skeleton.html`](./sessions/session3-lean-skeleton.html) | Mini-curso **S3** · 60 min · Lean 4 skeleton and AXLE Issue #12 |
| [`sims/helical-attractor.html`](./sims/helical-attractor.html) | Standalone Three.js (r128) 3D simulator for the dm³ flow |

---

## Where Chapter 10 sits in the series

*Principia Orthogona* is a five-volume series on **Generative Temporal Contact Theory (GTCT)**:

| G | Volume | Substrate | Output |
| --- | --- | --- | --- |
| G¹ | Vol I · GOMC | Abstract operator algebra | The orthogonal operator |
| G² | Vol II · TOGT | Contact geometry on 3-manifolds | Contact fixed point, g₃₃ = 33 |
| G³ | Vol III · Mini Beast | Biological instantiations | Living form |
| G⁴ | **Vol IV · GTCT-T1 (IMPA Edition) — this volume** | **Temporal contact** | **Helical attractors on the dm³ basin** |
| G⁵ | Vol V · AXLE | Formal verification in Lean 4 | Complete Completeness |

**Chapter 10 of Vol IV** is the self-contained exposition of the *dm³ operator* — a toy-model realisation of the full GTCT operator chain `G = U ∘ F ∘ K ∘ C` on the standard contact 3-manifold `(ℝ³, α)` with `α = dz − r² dθ` (polar coordinates). The helical attractors in the title are the orbits `Γ = {r = 1}` of the induced Reeb-like flow; the asymmetric basin described in Theorem 2.1 is the geometric object the entire course is organised around.

Parent publication (for libraries, retailers, metadata ingestion): *Principia Orthogona — IMPA Edition*, G6 LLC · Newark, NJ · 2026. Bowker-registered ISBNs identify the parent compilation, not this chapter:

| Format | ISBN-13 |
| --- | --- |
| Paperback | 979-8-9954416-7-0 |
| Hardback  | 979-8-9954416-8-7 |
| eBook (Adobe PDF) | 979-8-9954416-9-4 |

---

## The mini-curso — what each piece contains

Three 60-minute sessions, one embedded 3D simulator, one numerical lab, and a short primer chapter used as warm-up reading.

### Warm-up · [`chE-gtct.html`](./chE-gtct.html) — *GTCT for Everyone*

A standalone primer, originally written for Book 3 (*Mini-Beast*) and reused here. It answers "what is an axiom?" then states the **nine axioms of GTCT**, walks through the **twelve operators** as a map of the generative cycle (with an embedded *GTCT Time Machine — G Orbit Simulator*), states the **four main theorems**, reframes time through **five properties** (Order, Novelty, Rhythm, Irreversibility, Emergence), and closes with the **Seed Theorem** and a Prompt Panel. Participants read it before Session 1 and use it between sessions as a one-page reference card. This page is also the canonical home of the repo's *Principia Orthogona* dark-teal design language — every other page inherits from it.

### Session 1 · [`sessions/session1-contact-geometry.html`](./sessions/session1-contact-geometry.html) — *Contact Geometry and the dm³ System*

The first session builds the geometric stage. Outline:

1. **What is a contact 3-manifold?** The pair `(M, α)` with `dα|_{ker α}` non-degenerate.
2. **The prototype — ℝ³ with α = dz − r² dθ.** The polar Darboux model the whole chapter is written on.
3. **The Reeb vector field** and its flow on the prototype.
4. **The dm³ system.** The full three-equation ODE is exhibited and dissected:
   - *Anatomy of the radial equation* `ṙ = r(1 − r²) + ε(r − 1) e^(−z)` — logistic term plus coupling.
   - *Anatomy of the vertical equation* `ż = r² − ε(r − 1)² e^(−z)` — gain minus decay.
   - *Why ε matters* — the coupling strength that makes the basin asymmetric (fixed at ε = 2 throughout Vol IV).
5. **Live simulator** — the embedded `helical-attractor.html` runs inside the page with a guided 15-minute experiment block.
6. **Hand-off to Session S2.**

### Session 2 · [`sessions/session2-theorem-basin.html`](./sessions/session2-theorem-basin.html) — *Theorem 2.1 and the Asymmetric Basin*

The mathematical core of the mini-curso.

1. **Theorem 2.1 — Exponential Contraction to Γ.** For the dm³ system with ε = 2, if `|r(0) − 1| < 1/3` and `z(0) ≥ log 2`, then `r(t) → 1` exponentially with rate `μ = −2`; otherwise the orbit may escape and the basin is asymmetric.
2. **Linearisation at r = 1.** Where the `−2` comes from in one line — the derivative of the logistic term evaluated on Γ.
3. **The Grönwall sketch** (three steps):
   - Step 1 — the scalar inequality `d/dt |r − 1|² ≤ … `.
   - Step 2 — bound the bad pieces (the coupling term and the `e^(−z)` factor).
   - Step 3 — solve for `ε₀`, the allowed displacement.
4. **Table 1** — the numerical check at DOP853 precision (rtol = 10⁻¹⁰, atol = 10⁻¹²), with a `Read the table` segment and a "What goes wrong geometrically" explanation of inner-basin failure.
5. **Live numerical demo** — the embedded simulator is used to watch the inner edge of the basin in real time (8-minute guided demo).
6. **Hand-off to Session S3.**

### Session 3 · [`sessions/session3-lean-skeleton.html`](./sessions/session3-lean-skeleton.html) — *Lean 4 Skeleton and AXLE Issue #12*

Session 3 migrates the Session 2 proof into Lean 4 and leaves the participant holding a specific open conjecture.

1. **Why Lean 4, why now** — the formalisation case for contact dynamical systems.
2. **`AXLE/Chain.lean` — the operator chain.** Reading the Lean structure for `G = U ∘ F ∘ K ∘ C`.
3. **`SpiralReturn` — formalising Theorem 2.1.** What is already proved (the exponential-contraction *conclusion* and the linearisation *base case*) versus what is still `sorry` (the inductive step).
4. **AXLE Issue #12 — the Lipschitz estimate.** The specific `sorry` the participant leaves the session tracking; two already-closed issues are shown as precedent for how such estimates have been discharged before in AXLE.
5. **How to contribute.** A short list of Mathlib files to read before opening a PR, and a reference sim for grounding the proof in the flow it describes.
6. **What the course did.**

### Interactive simulator · [`sims/helical-attractor.html`](./sims/helical-attractor.html)

A standalone Three.js (r128) visualisation of the dm³ flow on `(ℝ³, α)`. Orbits are integrated in the browser; the helical attractor `Γ = {r = 1}` is rendered as a persistent reference curve; drag to rotate, scroll to zoom, drop initial points to watch them fall into the outer basin (contraction to Γ) or the inner basin (collapse). Embedded inside Sessions 1 and 2 as well; also reachable from the landing page's 3D hero.

### Numerical lab · [`labs/dm3_numeric.py`](./labs/dm3_numeric.py)

A single-file SciPy script that integrates the dm³ ODE with `solve_ivp(method='DOP853')` at rtol = 10⁻¹⁰ / atol = 10⁻¹², fits the empirical contraction rate `μ̂` as the slope of `log|r(t) − 1|` on `t ∈ [5, 20]`, records the early-time `ż` diagnostic, and prints the canonical Table 1.

```bash
pip install numpy scipy
python3 labs/dm3_numeric.py            # prints Table 1 to stdout
python3 labs/dm3_numeric.py --plot     # also saves dm3_trajectories.png (requires matplotlib)
```

Expected output: outer-basin rows converging to `μ̂ = −2`, inner-basin rows showing solver collapse when the integration falls off the asymmetric edge.

---

## AXLE — formal-verification companion

The Lean 4 engine that proves the dm³ operator's structural invariants lives in a sibling repository:

- **AXLE** (Lean 4 + Mathlib4): <https://github.com/TOTOGT/AXLE> — contains `Chain.lean` and `SpiralReturn.lean` referenced in Session 3.
- **Issue #12** — the open Lipschitz estimate: <https://github.com/TOTOGT/AXLE/issues/12>

Companion theory repository — **GTCT** (working notes and drafts): <https://github.com/TOTOGT/GTCT>.

---

## XII Bienal da SBM 2026 — submission artefacts

The same Chapter 10 material is submitted in four SBM categories. Files live in [`submissions/`](./submissions/). The SBM *Eixo Temático* mapping used here:

| Eixo | Theme |
| --- | --- |
| T1  | Belos Problemas e Belas Soluções |
| T2  | História da Matemática |
| T10 | Geometria e Topologia |

| Artefact | Category | Eixo |
| --- | --- | --- |
| `submissions/XII_BM_MINICURSO_T1_Pablo_Grossi_preview.pdf` | Minicurso (preview) | T1 |
| `submissions/PO_10_Pablo_Grossi.{tex,pdf,…}` | Pôster | T10 (Geometria e Topologia) |
| `submissions/PO_10_Pablo_Grossi_source.pdf` | Pôster typesetting source reference | — |
| `submissions/PO_10_T2_Pablo_Grossi.{tex,pdf,…}` | Pôster (second eixo) | T2 |
| `submissions/CO_T2_Pablo_Grossi.{tex,pdf,…}` | Comunicação Oral | T2 |
| `submissions/OF_T2_Pablo_Grossi.{tex,pdf,…}` | Oficina | T2 |

**SBM filename convention:** `{PREFIX}_T{N}_{FirstName}_{LastSurname}.pdf` where `PREFIX ∈ {MC, OF, PO, CO, Ex}` and `{N}` is the Eixo number — the `PO_10` / `MINICURSO` filenames will need the `T` prefix inserted before upload (e.g., `PO_T10_Pablo_Grossi.pdf`, `MC_T1_Pablo_Grossi.pdf`).

---

## Repository layout

```
3M/
├── README.md                                           this file
├── index.html                                          publisher landing page
├── chE-gtct.html                                       Chapter E · GTCT for Everyone (from Book 3)
│
├── assets/
│   ├── book-cover.png                                  rendered cover (from preprint p.1)
│   ├── book-cover-raw-01.png                           raw render
│   └── book-cover-sm.png                               thumbnail
│
├── sessions/                                           three 60-min mini-curso handouts
│   ├── session1-contact-geometry.html                  S1 · contact geometry + dm³
│   ├── session2-theorem-basin.html                     S2 · Theorem 2.1 + Table 1
│   └── session3-lean-skeleton.html                     S3 · Lean 4 / AXLE Issue #12
│
├── sims/
│   └── helical-attractor.html                          Three.js r128 3D dm³ simulator
│
├── labs/
│   └── dm3_numeric.py                                  DOP853 reproduction of Table 1
│
└── submissions/                                        XII Bienal SBM submission bundles
    ├── CO_T2_Pablo_Grossi.{tex,pdf,aux,log,out}        Comunicação Oral · T2
    ├── OF_T2_Pablo_Grossi.{tex,pdf,aux,log,out}        Oficina · T2
    ├── PO_10_Pablo_Grossi.{tex,pdf,aux,log,out}        Pôster · T10
    ├── PO_10_Pablo_Grossi_source.pdf                   typesetting source reference
    ├── PO_10_T2_Pablo_Grossi.{tex,pdf,aux,log}         Pôster · T2
    ├── XII_BM_MINICURSO_T1_Pablo_Grossi_preview.pdf    Minicurso preview · T1
    └── assets/                                         submission-local assets
```

---

## Design language

All pages inherit the *Principia Orthogona* dark-teal publisher aesthetic defined in `chE-gtct.html`.

- **Palette:** `--t #2dd4bf` (teal), `--bg #060f0e` (near-black), `--bg2 #0c1a18` (panel), `--gd #c9a84c` (accent gold).
- **Operator-chain colours:** `--c #4a9eff` (Contract), `--k #e05a3a` (Curvature), `--f #50c878` (Filter), `--u #c084fc` (Unfold).
- **Typography:** Georgia serif for prose, Courier New for code and chrome.
- **Math rendering:** MathJax 3 (`tex-svg`). 3D: Three.js r128. All from CDN — no build step, no install.

---

## License

© 2026 Pablo Nogueira Grossi — G6 LLC.
Code: MIT License. Text, figures, and session handouts: CC BY 4.0 unless otherwise noted.
ORCID: [0009-0000-6496-2186](https://orcid.org/0009-0000-6496-2186).

---

*Mathematics is a language. These theorems have been proved in every language simultaneously.*
*A matemática é uma língua. Os teoremas acima foram provados em todas as línguas simultaneamente.*

**C → K → F → U → ∞**
