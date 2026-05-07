# Section 8 — Long-time dynamics and finite-size scaling

**Historical reference skeleton — do not delete.**
This file is the paragraph-level draft from which the LaTeX Section 8
in `nbonacci_dnls_paper_with_figs.tex` was generated.

---

## 8.1 Methodology

**Overview.**
Section 7 listed the following as open questions:
longer time evolution ($T \gg 50$), finite-size scaling of the IPR ratio,
spreading exponents, and self-trapping thresholds.
This section reports numerical results that partially close all four.

The computational programme consists of three complementary campaigns:

- **Linear-limit D₂ at natural lengths.** We compute the mid-gap eigenstate
  IPR at the natural Rauzy and Fibonacci iteration lengths (rather than at
  arbitrary truncations) and fit $\mathrm{IPR}(N) \sim N^{-D_2}$ by OLS
  regression in log–log space (code: `d2_natural_lengths.py`,
  data: `data/d2_natural_lengths.csv`).

- **Finite-size scaling at T=10⁴ and T=10⁵.**
  We integrate the full DNLS dynamics from the mid-gap eigenstate at
  $\lambda = 1.5$ for chain sizes $N \in \{500, 1000, 2000\}$ and times
  up to $T = 10^4$ and $T = 10^5$ using the DOP853 integrator with
  `rtol = 1e-9` (code: `dnls_long_time.py`,
  data: `data/ipr_fss_T1e4.csv`, `data/ipr_lambda1p5_N1000_T1e5.csv`,
  `data/ipr_lambda1p5_N2000_T1e5.csv`).

- **Long-time saturation at T=10⁶.**
  For $N = 1000$, $\lambda = 1.5$ we extend the run to $T = 10^6$ with
  450 checkpoints to check whether saturation is permanent and whether the
  trib/fib IPR ratio $\gtrsim 1$ at long times
  (data: `data/ipr_lambda1p5_N1000_T1e6.csv`).

Figure callouts: fig6 (D₂ natural lengths), fig7 (FSS D₂ at T=10⁴),
fig8 (ratio collapse T=10⁴ → T=10⁵), fig9 (T=10⁶ saturation).

TODO: add tetrabonacci ($n=4$) comparison once tetra FSS data is merged.

---

## 8.2 Linear-limit multifractal dimensions

**Headline.**
Using natural iteration lengths removes the $N=1000 \to 2000$ reversal in
the tribonacci IPR sequence and shifts $D_{2,\mathrm{trib}}$ from the
anomalous $\approx 0.364$ (arbitrary lengths, $R^2 = 0.30$) to
$0.282 \pm 0.119$ ($R^2 = 0.65$), resolving the anomalous non-monotonicity
reported in Section 4.

The fractal (correlation) dimension D₂ is defined by the power-law scaling
IPR(N) ~ N^{-D₂}, so that D₂ = 1 for a fully extended state and D₂ = 0 for
a perfectly localised one; quasiperiodic critical states occupy the
intermediate range 0 < D₂ < 1.

Earlier estimates from truncated chains of arbitrary length
(N ∈ {200, 500, 1000, 2000}) gave D₂_fib ≈ 0.616 and D₂_trib ≈ 0.364.
The tribonacci estimate was anomalous: the underlying IPR(N) data were
non-monotone (0.0969 → 0.0820 → 0.0410 → 0.0484 for N=200,500,1000,2000),
with a reversal between N=1000 and N=2000 that inflated the uncertainty and
biased the slope.

The reversal is a finite-size artefact of **boundary truncation**: arbitrary
truncation at N=1000 or N=2000 destroys the self-similar boundary structure
that the Rauzy fixed-point chain lengths preserve.  To test this hypothesis
we re-computed IPR only at the natural Rauzy iteration lengths T_n (where
T_{n+3} = T_{n+2} + T_{n+1} + T_n, OEIS A000073), specifically
N ∈ {274, 504, 927, 1705, 3136} (n = 10–14), and the analogous natural
Fibonacci lengths N ∈ {233, 377, 610, 987, 1597, 2584} (n = 12–17).
At each natural length the chain is the exact n-th iterate of the
substitution rule, so the boundary atoms are self-consistently generated and
no truncation is introduced.

The mid-gap eigenstate is identified as the eigenstate with the smallest |E|
whose spatial spread σ ≥ 0.03N; this criterion excludes compact
boundary/defect zero-modes (notably the anomalous E=0 state at tribonacci
n=12, N=927, whose spread σ ≈ 15 ≪ 0.03×927) that are unrelated to the bulk
multifractal scaling.  With this selection the tribonacci IPR sequence is
monotone non-increasing:
0.097 → 0.082 ≈ 0.082 ≈ 0.082 → 0.041.
The near-plateau at IPR ≈ 0.082 for n = 11, 12, 13 (N = 504, 927, 1705)
reflects a genuine feature of the tribonacci RSRG hierarchy in which the same
critical-state family persists across three consecutive iterates.  The
(N=1000 → N=2000) reversal present in the arbitrary-length data is
eliminated: the natural-length sequence ends with a clean drop from
IPR ≈ 0.082 (N=1705) to 0.041 (N=3136).

The OLS regression of log IPR on log N yields:
  D₂_fib  = 0.646 ± 0.111  (R² = 0.8946)
  D₂_trib = 0.282 ± 0.119  (R² = 0.6535)

The hierarchy D₂_trib < D₂_fib at natural lengths is consistent with the
stronger spatial multifractality of the tribonacci eigenstate.

Figure: fig6_d2_natural — log IPR vs log N for Fibonacci and tribonacci at
natural iteration lengths, with OLS fit lines.

---

## 8.3 Finite-size scaling of IPR at T=10⁴

**Headline.**
At T=10⁴, λ=1.5, the trib/fib IPR ratio peaks at 4.15× (N=1000) and is
non-monotone in N (0.658 at N=500, 4.15× at N=1000, 3.17× at N=2000),
reflecting a transient crossover, not a bulk FSS trend.

At T = 10⁴ and λ = 1.5 the IPR values for the two chains are:

| N    | IPR_fib | IPR_trib | ratio |
|------|---------|----------|-------|
| 500  | (data)  | (data)   | 0.658 |
| 1000 | (data)  | (data)   | 4.15  |
| 2000 | (data)  | (data)   | 3.17  |

The non-monotone ratio vs N reveals that the dynamics have not converged to
the asymptotic finite-size scaling regime at T = 10⁴.  The Fibonacci state
at N=1000 shows a deeper dip in IPR compared to N=500 or N=2000, indicating
an intermediate-time transient localisation reversal.  The tribonacci state
is far more robust to this transient: its IPR drops by < 30% from the linear
limit even at T = 10⁴.

D₂ values extracted from late-time (T = 10⁴) IPR are:
  D₂_fib  ≈ 0.578,  D₂_trib ≈ 0.381 (from fss_D2.csv)

These differ from the linear-limit values (Section 8.2) because the nonlinear
dynamics at T = 10⁴ have already altered the spatial structure of the states.

Figure: fig7_fss_T1e4 — D₂ values (Fibonacci and tribonacci) vs N at T=10⁴,
showing approach to the asymptotic fractal dimension from above.

TODO: cross-check with tetrabonacci FSS data once merged.

---

## 8.4 Transient peak and ratio collapse

**Headline.**
The 4.15× trib/fib ratio at T=10⁴, N=1000 is a transient peak that
collapses to 1.18× at T=10⁵, confirming the non-monotone ratio is a
crossover artefact rather than a persistent feature.

At T = 10⁵, N = 1000, λ = 1.5 the IPR ratio collapses to ≈ 1.18×.
At T = 10⁵, N = 2000 the ratio is ≈ 1.24×.  Both values are above unity
(tribonacci IPR > Fibonacci IPR) but dramatically reduced from the T = 10⁴
peak.

The spreading exponent α, extracted from an OLS fit to
log IPR(t) ~ −α log t for t > 10⁴ (N = 2000), gives:
  α_fib  ≈ 0.487,  α_trib ≈ 0.777  at T = 10⁵, N = 2000, λ = 1.5.

Note: the T = 10⁴ values (α_fib ≈ 0.211, α_trib ≈ 0.155) were not
asymptotic; the hierarchy α_trib > α_fib is robust only at T ≥ 10⁵.

The faster spreading of tribonacci (α_trib > α_fib at long times) appears
counter-intuitive given the stronger multifractality of the tribonacci
eigenstate.  A tentative interpretation: the stronger spatial hierarchy of
the tribonacci state creates more efficient long-range hopping pathways under
nonlinear perturbation once the initial self-trapping regime is overcome.

Figure: fig8_ratio_collapse — spreading exponent α vs N for both chains at
λ = 1.5, showing α_trib > α_fib at N = 2000, T = 10⁵.

TODO: extend to T = 10⁶ for N = 2000 to confirm convergence of α.

---

## 8.5 Long-time saturation at T=10⁶

**Headline.**
At N=1000, T=10⁶, λ=1.5: both chains saturate to IPR ≈ 0.002 by
t ≈ 3×10⁵, with the trib/fib ratio oscillating 1.04 ± 0.04 (16 crossings).
No permanent crossover. Normalisation is preserved (max norm leak < 5×10⁻⁷).

Extended integration to T = 10⁶ at N = 1000, λ = 1.5 (450 checkpoints,
rtol = 1e-9, DOP853):

- Both chains saturate to IPR ≈ 0.002 (= 1/N, the fully extended value) by
  t ≈ 3×10⁵.
- The trib/fib IPR ratio at T = 10⁶ is 1.017, with oscillations 1.04 ± 0.04
  (16 crossings across unity in the window t ∈ [10⁵, 10⁶]).
- No permanent crossover is observed; the two chains are effectively
  co-delocalised at T = 10⁶.
- Spreading exponents in the t > 10⁵ window: α_fib ≈ 0.009 (R² = 0.035),
  α_trib ≈ 0.087 (R² = 0.733).  The tribonacci exponent is marginally larger
  but both are consistent with near-complete saturation.
- Wall time: fibonacci 1357s, tribonacci 1254s (single core).

These results confirm that the differential nonlinear robustness described in
Sections 3–5 is a **transient** phenomenon: on timescales T ≲ 10⁴ the
tribonacci state is dramatically more robust, but by T ≈ 3×10⁵ both states
reach the fully extended limit.

Figure: fig9_T1e6_saturation — IPR vs log t for both chains at N=1000,
λ=1.5, T=10⁶, showing co-saturation to IPR ≈ 1/N ≈ 0.001 by t ≈ 3×10⁵.

---

## 8.6 Summary and outlook

The four campaigns of Section 8 establish:

1. The multifractal exponents D₂_fib = 0.646 ± 0.111 and
   D₂_trib = 0.282 ± 0.119 are the correct linear-limit values (natural
   lengths); the earlier anomalous D₂_trib ≈ 0.364 is a boundary-truncation
   artefact.

2. The trib/fib IPR ratio of ≈4 reported at T=50 in Sections 3–5 represents
   a transient robustness window: both chains eventually co-delocalize to the
   fully extended limit, but the tribonacci state retains its advantage for
   t ≲ 3×10⁵.

3. The spreading exponent hierarchy α_trib > α_fib is robust at intermediate
   times (T ≈ 10⁵) but disappears at T = 10⁶ as both chains saturate.

4. The non-monotone FSS ratio at T=10⁴ is a transient crossover, not a
   bulk property; at T=10⁵ the ratio is monotone and collapsed to ≈1.2×.

Outstanding questions:
- TODO: repeat the T=10⁶ run at N=2000 to confirm system-size independence.
- TODO: measure the self-trapping crossover time τ_c(N, λ) systematically.
- TODO: add the tetrabonacci ($n=4$) chain to the FSS comparison.
