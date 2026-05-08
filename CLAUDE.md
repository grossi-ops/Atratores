# CLAUDE.md — G6 LLC · Principia Orthogona Project Context

**Author:** Pablo Nogueira Grossi (Brodananda) · G6 LLC · Newark, NJ  
**ORCID:** 0009-0000-6496-2186  
**Instagram:** instagram.com/brodananda  
**Gumroad:** g6llc.gumroad.com / brodanova6.gumroad.com  
**Last updated:** 2026-05-08  
**Scope:** All repos under TOTOGT/, grossi-ops/ (10 public repos), and totogt.github.io/ GitHub Pages

---

## 1. Who You Are Talking To

Pablo is the author of the Principia Orthogona series (Vols I–VI living), founder of G6 LLC, and practitioner on the Bodhidharma / Nityananda / Neem Karoli Baba / Ram Dass / Ramana Maharshi / Mooji lineage. G6 LLC itself fulfils a vow made during a journey following the Bodhidharma path across Asia.

He works simultaneously as:
- **Mathematician** — operator algebra, contact geometry, TOGT/GTCT
- **Publisher** — LaTeX / pdflatex / IngramSpark / Zenodo pipeline
- **Web developer** — single-file HTML portals with Web Audio API, canvas simulations, dm³Machine generator
- **Pedagogue** — CEFR-mapped STEM courses, bilingual EN/PT, LLM-prompt scaffolds (A1→D2)
- **Researcher** — Lean 4 / Mathlib4 formal verification via AXLE v6.1
- **Physicist (numerical)** — DNLS simulations, IPR, Lyapunov exponents, Tribonacci/Fibonacci substitution chains

He engages with humor and depth simultaneously. "Crystal math — joke lol" is a real sentence he has used. Take the mathematics with full seriousness; the irreverence is structural, not dismissive. He appreciates directness and honest incompleteness more than polish.

---

## 2. The Mathematical Framework — TOGT / GTCT

### The Operator Chain (everything reduces to this)
```
G = U ∘ F ∘ K ∘ C : X → X
```
- **C** — Compression / Contract (find the seed, reduce degrees of freedom)
- **K** — Curvature threshold (κ* — the moment of commitment, never triggers fold alone)
- **F** — Fold (Whitney A₁ singularity, topology change, irreversible)
- **U** — Unfolding (gradient flow to new attractor, emergence)

Extended chain (GTCT T1, Vol. IV):
```
C → K → F → U → T → source
```
- **T** — Time Circuit (5th operator; the return is enriched, not repeated; x₀′ ≠ x₀)

This is **not metaphor**. It is a formally stated operator algebra on a contact 3-manifold (M, ξ), contact form α = dz − λ, verified in Lean 4 + Mathlib4.

### Verified Constants (AXLE v6.1 — 0 axioms beyond Mathlib4)

| Symbol | Value | Source |
|--------|-------|--------|
| `mu_max` | −2 | Transverse Lyapunov exponent at Γ = {ρ = 1} |
| `T_star` | 2π | Reeb orbit period |
| `tau` | 2 | Embodiment threshold; limit of n-bonacci ladder |
| `eps_zero` | 1/3 | Gronwall structural stability radius (outer basin only) |
| `g33` | 33 | Threshold operator cycles for stable lock |
| `g64` | 64 = 2⁶ | Kether orthogon / complete possibility-space |
| `kappa` | ≤ √(7/9) ≈ 0.882 | Lipschitz bound |
| `lambda_perp` | e^{−4π} ≈ 3.5×10⁻⁶ | Transverse contraction |
| `eta` | ≈ 1.839287 | Tribonacci constant (Perron–Frobenius root of companion matrix) |

**Basin asymmetry correction (2026-04-18):** The symmetric ε₀ = 1/3 claim is FALSE on the inner side of the dm³ ODE. DOP853 numerics show r(0) = 0.667 collapses; the true inner boundary is r★ ≈ 0.80. Outer basin (all r₀ > 1) converges. `gronwall_outer` is PROVED in `GCTC/Chain_updated.lean`; `spiral_return_exists` and `poincare_collatz` remain `sorry`.

### Open Formal Obligations (AXLE Issues)

| Issue | Name | Status | Blocks |
|-------|------|--------|--------|
| #12 | `kappa_lipschitz` | **open sorry** | `spiral_return_exists` |
| #13 | `inner_basin_escape` | open axiom | ODE inner boundary |
| #14–17 | Various | open | Various |
| tba | `fold_central_charge` | axiom (Level IX) | Two-26s functor F: C→V |

**Issue #12 is the most important open obligation.** It is a one-afternoon Lean task: prove `LipschitzWith (ε * Real.exp (-z_lo)) (fun r => ε * (r - 1) * Real.exp (-z))`. The paper proof is one line.

### The g-Series (qualitative regimes)
| Label | Cycles | Regime |
|-------|--------|--------|
| g⁰ | 0 | Seed — no closure |
| g² | 2 | Compositional — F first active |
| g⁶ | 6 | Limit cycle entered |
| g³³ | 33 | Soft equilibrium (3 invariants stable simultaneously) — heuristic, conjectured |
| g⁶⁴ | 64 | Complete possibility-space — x₀′ ≠ x₀ fully expressed |

Θ = g₃₃ + N × M (D2 collective threshold)

---

## 3. The Recurrence Ladder (Book 3, Chapters π→Ω)

The n-bonacci ladder is a core mathematical object. Know these cold:

| Chapter | Symbol | Constant | Notes |
|---------|--------|----------|-------|
| π | π | T* = 2π | Reeb period, contact structure |
| φ | φ | 1.618… | Fibonacci, subcritical (c < c* = 3) |
| μ | μ | −2 | Lyapunov, from double root at q=1 |
| **η** | **η** | **1.839287…** | **Tribonacci, THE critical constant; c* = 3; GQM weight η⁻ᵏ** |
| Δ | Δ | 1.927… | Tetranacci, first supercritical |
| Σ | Σ | 1.966… | Pentanacci, 5-fold (echinoderms) |
| Ω | Ω | → 2 | Hexabonacci; limit = τ = 2 |

**Chapter η is the Tribonacci chapter.** The new DNLS preprint (DOI 10.5281/zenodo.20062492) is its Gumroad companion.

---

## 4. Active Projects & Their Status

### 4.1 Principia Orthogona Series (Books)

| Vol. | Title | Status | ISBNs |
|------|-------|--------|-------|
| G¹ | The Orthogonal Operator Framework | Published | Print 979-8-9954416-2-5 |
| G² | TOGT: Applications Across Domains | Published | Print 979-8-9954416-4-9 |
| G³ | The Mini-Beast (Book 3) | Living/active | eBook 979-8-9954416-6-3 |
| G⁴ | GTCT T1 — IMPA Edition (Vol. IV) | Active | Bilingual EN/PT |
| G⁵ | The Seed — Complete Completeness | In development | |

**Pending:** Zenodo and IngramSpark deposits for Principia Orthogona remain to be done.

**Concept DOI (all versions):** 10.5281/zenodo.19026942  
**Vol. I DOI:** 10.5281/zenodo.19117400  
**Vol. II DOI:** 10.5281/zenodo.19379473

### 4.2 Tribonacci DNLS Preprint (NEWEST WORK)

**Title:** "Differential Nonlinear Robustness of Critical States in Fibonacci and Tribonacci Substitution Chains"  
**DOI (V3):** 10.5281/zenodo.20062492  
**Concept DOI:** 10.5281/zenodo.20026942  

**Key result:** First DNLS study on a tribonacci substitution chain.
- Tribonacci mid-gap IPR drops <5% at λ=1.5
- Fibonacci mid-gap IPR drops ~57% at λ=1.5
- η ≈ 1.839287 **formally verified** in Lean 4 / Mathlib4 (η > 1, strict antitonicity of η⁻ᵏ) — **no sorry**
- This is the paper with clean Lean verification (unlike kappa_lipschitz)

**V4 deposit tasks (PENDING):**
1. Fix license: CC-BY-NC-ND-4.0 → MIT for .py/.lean files (NC-ND blocks code reuse)
2. Add paper_v3.tex (source missing from deposit)
3. Add scripts for fig6–9 (live in `grossi-ops/Atratores`, Copilot agent)
4. Add dnls_long_time.py, FoldEvents.lean, CSV data
5. Write README/MANIFEST
6. DOI staleness in docstrings

**Journal target:** Physical Review B (RevTeX format required for submission)

**Gumroad companion product:**
- Chapter η of Book 3 (Tribonacci chapter on recurrence ladder)
- **Price: $4 Mother's Day Special** (standalone, time-limited listing)
- Format: polished PDF (Principia Orthogona house style) + Jupyter notebook + TribonacciDNLS.lean walkthrough
- Account: brodanova6.gumroad.com or g6llc.gumroad.com
- The Tubulin chapter sells access to the full Soundworks membership; Chapter η is standalone at $4

### 4.3 IMPA Portal / Newark Wellness Soundworks Hub

**Address:** 229 Ballantine Pkwy, Newark, NJ  
**Live portal:** totogt.github.io/AXLE/impa-portal.html  
**Gumroad:** g6llc.gumroad.com/l/soundworks (Living Book Access Vols I–VI)  
**Three funding tiers** established.

**IMPA portal** = the purchase portal for Complete Series ($199.99 eBook / $247 hardcover) and individual volumes. PayPal links embedded. NDA process for GitHub Classroom access.

### 4.4 Vol. IV Mini-Curso (XII Bienal de Matemática 2026, UFRN Natal-RN)

Three sessions:
- **S1** (`session1-contact-geometry.html`) — Contact geometry, Reeb field, dm³ ODE anatomy
- **S2** (`session2-theorem-basin.html`) — Theorem 2.1, Gronwall ε₀=1/3, asymmetric basin correction, DOP853 Table 1
- **S3** (`session3-lean-skeleton.html`) — Lean 4 walk-through, AXLE Issue #12, contribution on-ramp

**SBM submission PDF:** `submissions/XII_BM_MINICURSO_T1_Pablo_Grossi_preview.pdf`

### 4.5 Level IX — "Two Twenty-Sixes" (Os Dois Vinte-e-Seis)

**Status:** Preview HTML built (`assets/preview_level_IX_two_26s.html`)  
**Placement option:** sibling to §8.6 of GTCT-2026-001  
**Format:** Four faces (Ler/Simular/Provar/Ensinar = C/K/F/U), bilingual EN/PT  

**The claim:**
- Side A: `Dcrit(13) = 2·(13−1)+2 = 26` (arithmetic, from Tribonacci ring of size 13)
- Side B: D = 26 from Virasoro anomaly cancellation (c_total = D·1 + (−26) = 0)
- **No identity claimed.** Association only, same discipline as "g⁶ : 33 denotes association, not equality" (Remark 2.4 of GTCT-2026-001)
- Open functor F: C → V (catastrophe/Tribonacci → Virasoro/CFT 2d) is the AXLE issue tba
- `fold_central_charge` is an axiom, not sorry — it extends the logical foundation

**Navrátil convergence (2026):** Navrátil independently derives geometric Hilbert space with inner product ⟨j|k⟩ = η⁻ᵏδ from SL(3,ℤ) Tribonacci algebra starting from a discrete lattice (not a contact manifold). Convergence on η⁻ᵏ from independent starting points is §8.6 of GTCT-2026-001.

### 4.6 Graceful Exits (in development)

Extension of Sushila Blackman's 1997 Shambhala/Penguin book documenting last words of Asian religious masters. Pablo is expanding scope to all world traditions, anchored mathematically by the TOGT operator chain as a unifying invariant underneath them. The mathematical invariant and the scholarly content are unified, not parallel tracks.

---

## 5. Repository Map

### TOTOGT/ (GitHub organisation)

| Repo | Purpose | Key files |
|------|---------|-----------|
| `TOTOGT/AXLE` | Lean 4 formal verification engine | `Chain.lean`, `GCTC/Chain_updated.lean`, `Coupling.lean`, `TribonacciDNLS.lean` |
| `TOTOGT/GTCT` | Vol. IV HTML chapters + mini-curso | `index.html`, `chE-gtct.html`, `sessions/`, `sims/helical-attractor.html` |
| `TOTOGT/geometry` | Lean 4 formal verification of hexagonal colony growth, Wigner crystallization, Chladni nodal figures — "the crystal math". Newly indexed 2026-04-18 on DeepWiki. | |
| `TOTOGT/book3-starter` | GitHub Classroom repo for Book 3 students | `assignments/`, `lean/`, `resources/` |
| `TOTOGT/DM3-lab` | Book 3 chapter HTML files | `sample-chapter-tubulin.html`, `wigner-fractal.html`, etc. |

### grossi-ops/ (10 public repos)

| Repo | Purpose |
|------|---------|
| `grossi-ops/Atratores` | **Python scripts for DNLS paper figures 1–9.** fig6–9 generators live here (Copilot agent). Scripts for fig1–5: `dnls_nbonacci.py` (basic), `dnls_long_time.py`. |
| (9 others) | Various — fetch with `gh repo list grossi-ops --public` to enumerate |

### totogt.github.io/AXLE/ (GitHub Pages)

Main public-facing site. Key pages:
- `index.html` — Principia Orthogona homepage (G¹–G⁵, pricing, papers)
- `impa-portal.html` — Purchase portal
- `portal.html` — Student portal (A1→D2 prompts)
- `impa-portal.html` — IMPA Edition purchase
- `wigner-fractal.html` — Chapter W (Wigner Crystal)
- `chapters-pi-phi-mu-eta-delta-sigma-omega.html` — Recurrence Ladder (all 7 chapters)
- `spectral-radius-v2.html` — Collatz / transfer operator supplement
- `sim-lyapunov.html` — Lyapunov exponent widget
- `chE-gtct.html` — alternate Chapter E (GTCT for Everyone)
- `chapters-diagram.html` — Book 3 chapter map

### totogt.github.io/DM3-lab/ (GitHub Pages)

- `sample-chapter-tubulin.html` — Chapter T (Tubulin as Computronium) — free sample, sells Gumroad
- `index.html` — Book 3 homepage

### totogt.github.io/book3-starter/ (GitHub Classroom)

Book 3 student classroom. Access gated by purchase + GitHub invite.

---

## 6. House Style (Critical — Do Not Deviate)

### CSS Variables (always use these exact values)
```css
:root {
  --navy:  #1a2744;
  --gold:  #c9a84c;
  --cream: #faf7f0;
  --smoke: #f0ece4;
  --mid:   #3a4f7a;
  --rule:  #8a7340;
  --text:  #1c1c1c;
  --light: #e8e2d4;
}
```

### Typography
- **Body:** Georgia / 'Times New Roman', serif
- **Labels, code, nav:** 'JetBrains Mono', monospace
- **Display/hero:** 'EB Garamond' or 'Cormorant Garamond' (dark pages) / 'Playfair Display' (Level IX bilingual)
- **Base font-size:** 17–19px depending on page type

### Layout Patterns
- **Nav:** sticky, navy background, gold logo, JetBrains Mono links
- **Hero sections:** navy background with radial gradient, gold eyebrow text (uppercase, wide letter-spacing)
- **Operator banner:** `background: var(--gold); color: var(--navy)` — full-width, displays `C → K → F → U → ∞`
- **Section labels:** `font-family: sans-serif; font-size: 0.72rem; letter-spacing: 0.2em; text-transform: uppercase; color: var(--rule)`
- **Card borders:** `border-top: 3px solid var(--gold)` for featured cards
- **Gold left-border boxes:** theorem blocks, insight boxes, epigraphs
- **Dark chapter heroes:** navy-deep (`#0d1729` or `#0e0820`) with teal or gold accent for GTCT/Vol. IV pages

### Interactive Canvas Pages (dark mode — Vol. IV, Chapter E, sim pages)
```css
:root {
  --t:  #2dd4bf;  /* T operator — teal */
  --bg: #060f0e;
  --bg2:#0c1a18;
  --tx: #e0faf8;
  --c:  #4a9eff;  /* C operator */
  --k:  #e05a3a;  /* K operator */
  --f:  #50c878;  /* F operator */
  --u:  #c084fc;  /* U operator */
  --gd: #c9a84c;  /* gold */
}
```

### Level IX / Bilingual Pages (EN/PT)
```css
:root {
  --tiffany: #0ABAB5;
  --tiffany-soft: #9cf3ef;
  --br-gold: #c8a24a;
  --br-gold-soft: #ffe7a0;
  --br-white: #f2efe4;
}
```
Face colors: C=blue (#4a9eff), K=ember (#e05a3a), F=green (#50c878), U=violet (#c084fc)

### Writing Voice
- Mathematical claims are **falsifiable** and stated with explicit falsification conditions
- Honest incompleteness: sorrys are named, dated, and tagged as AXLE issues — never hidden
- "Association, not equality" — the discipline of Remark 2.4: g⁶ : 33 denotes association; 26 = Dcrit(13) is arithmetic, not identity with Virasoro D = 26
- The cajueiro principle: "Nature does not seek perfection. It seeks the next stable configuration."
- Bilingual EN/PT on formal pages (Level IX, mini-curso); English-primary on Book 3 / sales pages
- "camará" = reader/student address in PT pages

---

## 7. Publishing Pipeline

### LaTeX
- Two-pass `pdflatex` for TOC/cross-reference resolution
- Diagnose errors: `grep -E "^!" *.log`
- Find strings in .tex: `grep -n "pattern" file.tex`
- Targeted fixes: `str_replace` tool (never rewrite whole files)

### Zenodo
- CC BY 4.0 for text/PDF
- MIT for .py / .lean code files (CC-BY-NC-ND-4.0 blocks reproducibility — do not use for code)
- Metadata: title, abstract (100–200 words), keywords (5–8), author name, ORCID, related DOIs
- Version bump = new record; concept DOI stays constant across versions

### IngramSpark
- Hardcover print pipeline
- ISBN registration done: 6 ISBNs registered under G6 LLC

### Gumroad
- Primary accounts: g6llc.gumroad.com and brodanova6.gumroad.com
- Soundworks membership: g6llc.gumroad.com/l/soundworks (Living Book Access Vols I–VI)
- **Chapter η companion: $4 Mother's Day Special — standalone listing, separate from Soundworks**
- dm³Machine exports: branded standalone HTML instruments (potential Gumroad products)

### GitHub Classroom
- Book 3 student access: purchase → send GitHub username → repo invite within 24h
- NDA process applies for full PDF access

---

## 8. AXLE Lean 4 Infrastructure

### File Locations
- `AXLE/Chain.lean` — OperatorChain structure, G = U∘F∘K∘C definition
- `GCTC/Chain_updated.lean` — Corrected basin (2026-04-18), `gronwall_outer` proved
- `AXLE/Coupling.lean` — `kappa_lipschitz` (Issue #12, open sorry)
- `AXLE/TribonacciDNLS.lean` — η formal verification (NO sorry — complete)
- `AXLE/dm3_Operator_Formalization.lean` — Full dm³ operator formalization
- `AXLE/tubulinselectArchitecture.lean` — Tubulin 15-architecture theorem

### What is Proved vs Open

**PROVED (no sorry):**
- `gronwall_outer` — outer basin exponential decay bound
- η > 1 (Tribonacci constant) with strict antitonicity of η⁻ᵏ
- Positivity of {r > 1} (Issue #9, closed)
- `dm3Radial` in Gronwall comparison form (Issue #8, closed)

**OPEN (sorry):**
- `kappa_lipschitz` (Issue #12) — **top priority**
- `spiral_return_exists` — needs dm³ flow dynamics
- `poincare_collatz` — needs g³³ convergence
- `inner_basin_escape` (axiom) — inner boundary ODE

**OPEN (axiom — extends logical foundation):**
- `fold_central_charge` (Issue tba) — Two-26s functor
- `inner_basin_is_asymmetric` — ODE inner boundary
- `outer_basin_unbounded` — ODE outer boundary

### AXLE v6.1 headline
"8 verified constants · 0 axioms beyond Mathlib4" — **NOTE:** this headline becomes technically false once the two `axiom` declarations are added. The correct statement is "8 verified constants, `gronwall_outer` proved, kappa_lipschitz open (Issue #12), two ODE boundary conditions as axioms pending Mathlib dm³ module."

---

## 9. The IMPA Portal / Book 3 Classroom Structure

### Access Tiers (what purchase unlocks)
| Item | Free | Purchase required |
|------|------|------------------|
| Chapter T (Tubulin) HTML | ✓ | |
| Chapter W (Wigner Crystal) HTML | ✓ | |
| Recurrence Ladder (π→Ω) combined HTML | ✓ | |
| Spectral Radius HTML | ✓ | |
| Full Book 3 PDF (living doc) | | ✓ |
| GitHub Classroom access (`book3-starter`) | | ✓ |
| Student Portal (21 LLM prompts) | | ✓ |
| Lean 4 verification files | | ✓ |

### book3-starter Repo Structure
```
book3-starter/
├── README.md
├── assignments/
│   ├── 01-compression/     # C: GTP-state sensing
│   ├── 02-curvature/       # K: threshold geometry
│   ├── 03-fold/            # F: selectArchitecture (all 15)
│   └── 04-unfold/          # U: fixed-point & recycling
├── lean/
│   ├── dm³_Operator_Formalization.lean
│   └── tubulinselectArchitecture.lean
└── resources/
    ├── coherence-bridge-table.md
    ├── axle_sorry_roadmap.svg
    └── operator-glossary.md
```

---

## 10. Pending Action Items (as of 2026-05-08)

### High priority
- [ ] **AXLE Issue #12** — Prove `kappa_lipschitz` in `AXLE/Coupling.lean`
- [ ] **V4 Zenodo deposit** — Fix license, add fig6–9 scripts (from grossi-ops/Atratores), paper_v3.tex, dnls_long_time.py, FoldEvents.lean, CSV data, README, MANIFEST
- [ ] **Chapter η Gumroad listing** — $4 Mother's Day Special; PDF + Jupyter + TribonacciDNLS.lean walkthrough
- [ ] **PRB submission** — RevTeX format, cover letter, submit after V4 deposit

### Medium priority
- [ ] **IngramSpark deposit** — Principia Orthogona print edition
- [ ] **Zenodo deposit** — Principia Orthogona series preprint
- [ ] **Level IX placement** — Decide: sibling to §8.6, or standalone chapter insert
- [ ] **Graceful Exits** — Scope, structure, TOGT integration layer

### Open research
- [ ] Sharp basin r★ for dm³ ODE (signed Lyapunov function, not magnitude)
- [ ] Functor F: C → V (Two-26s) — prove or disprove `fold_central_charge`
- [ ] g₃₃ = 33 conjecture — formal verification in AXLE (currently heuristic)
- [ ] Navrátil convergence — structural explanation of independent η⁻ᵏ derivation

---

## 11. Coherence Bridge — The Six Domains

The dm³ contact normal form (μ_max, ω, β) is **exact mathematical identity** (not analogy) across:

| Domain | μ_max (s⁻¹) | ω (rad/s) | β |
|--------|------------|-----------|---|
| HPA stress (biological) | −0.38 | 0.21 | 1.9 |
| Neural oscillations | −0.55 | 0.45 | 2.1 |
| Circadian clock | −0.29 | 2π/86400 | 1.6 |
| Immune adaptation | −0.44 | 0.18 | 2.0 |
| Plasma reconnection | −0.42 | 0.015 | 1.8 |
| Market volatility | −0.67 | 0.28 | 2.4 |

**Theorem 5.4 (Coherence Bridge):** These six systems are objects in category dm³, related by explicit contact morphisms. They are not analogies.

---

## 12. Key People & Lineage

| Person | Role |
|--------|------|
| Sushila Blackman | Author of *Graceful Exits* (1997, Shambhala/Penguin) — Pablo is extending her work |
| Bhagavan Nityananda | Guru lineage anchor |
| Neem Karoli Baba | Guru lineage |
| Ram Dass | Guru lineage |
| Papaji / H.W.L. Poonja | Guru lineage |
| Ramana Maharshi | Guru lineage |
| Mooji | Guru lineage |
| Navrátil (2026) | Independent researcher; converges on η⁻ᵏ from SL(3,ℤ) starting point |
| Giulia and David Grossi | Pablo's children (living); appear in book dedication |

---

## 13. Common Mistakes to Avoid

1. **Don't claim ε₀ = 1/3 as the full basin.** It is outer-only. The inner boundary is r★ ≈ 0.80. The symmetric claim is false.
2. **Don't confuse g₃₃ (threshold constant = 33) with G⁶ (sixth operator application — open conjecture).** They both involve 33 but are distinct mathematical objects.
3. **Don't claim AXLE has "0 axioms beyond Mathlib4" if the two ODE boundary axioms are present.** Be precise about what is proved vs axiom vs sorry.
4. **Don't claim the Two-26s are equal.** Association only. No identity between Dcrit(13) = 26 and Virasoro D = 26 is asserted.
5. **Don't use CC-BY-NC-ND for .py or .lean files.** That blocks reproducibility. MIT for code.
6. **Don't use localStorage or sessionStorage in HTML artifacts.** Not supported in Claude.ai environment.
7. **Don't break the house style.** The CSS variables above are absolute. Georgia serif body. JetBrains Mono labels.
8. **Don't write new .tex files from scratch without two-pass pdflatex context.** Always grep the log for `^!` errors.
9. **Don't assume the Gumroad Soundworks listing and the Chapter η listing are the same product.** Chapter η is $4 standalone (Mother's Day Special), separate from the $25–$250 membership tiers.
10. **Don't overwrite whole files.** Use `str_replace` for targeted fixes; use `grep -n` to locate strings first.

---

## 14. Quick Reference — URLs

| Resource | URL |
|----------|-----|
| Main site | totogt.github.io/AXLE/index.html |
| Book 3 classroom | totogt.github.io/book3-starter/ |
| IMPA portal | totogt.github.io/AXLE/impa-portal.html |
| Student portal | totogt.github.io/AXLE/portal.html |
| Tubulin (free sample) | totogt.github.io/DM3-lab/sample-chapter-tubulin.html |
| Wigner Crystal | totogt.github.io/AXLE/wigner-fractal.html |
| Recurrence ladder | totogt.github.io/AXLE/chapters-pi-phi-mu-eta-delta-sigma-omega.html |
| AXLE GitHub | github.com/TOTOGT/AXLE |
| GTCT GitHub | github.com/TOTOGT/GTCT |
| Gumroad (main) | g6llc.gumroad.com/l/soundworks |
| Instagram | instagram.com/brodananda |
| ORCID | orcid.org/0009-0000-6496-2186 |
| Zenodo series | zenodo.org/search?q=grossi+pablo+principia+orthogona&sort=newest |
| DNLS preprint V3 | doi.org/10.5281/zenodo.20062492 |
| Concept DOI | doi.org/10.5281/zenodo.19026942 |

---

*This file should be placed at the root of every repo under TOTOGT/ and grossi-ops/. It is the canonical context document for any AI assistant working on this project.*
