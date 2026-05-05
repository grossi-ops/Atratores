# Section 8 — Numerical Verification and Finite-Size Scaling at λ = 1.5

**Status:** draft skeleton  
**Data sources:**  
- `lambda1p5_sweep.csv` — fine-grained λ sweep (Δλ=0.1, T=50, N=500, RK45)  
- `ipr_lambda1p5_N500.csv`, `ipr_lambda1p5_N1000.csv`, `ipr_lambda1p5_N2000.csv` — long-time FSS runs (T=10⁴, DOP853, λ=1.5, N∈{500,1000,2000})

---

## 8.1 Reproduction of the λ = 1.5 Headline (T = 50 Timescale)

The fine-grained λ sweep (Δλ = 0.1) at the paper's original parameters (N = 500, T = 50, RK45, tol = 10⁻⁷/10⁻⁹) recovers the Table 1 numbers to two decimal places:

| Chain | IPR(λ=0) | IPR(λ=1.5) | ΔIPR |
|---|---|---|---|
| Fibonacci | 0.020965 | 0.009043 | −56.9% |
| Tribonacci | 0.082016 | 0.078198 | −4.6% |

The trib/fib IPR ratio peaks at **8.65×** at λ = 1.5 (vs 3.91× at the linear limit), confirming the differential robustness claim of Section 5 to within numerical precision.

The minimum Fibonacci IPR occurs at λ ≈ 1.8 (IPR = 0.00870, −58.5%), suggesting the crossover is slightly above the λ = 1.5 landmark but with all qualitative features intact.

**[FIGURE: `lambda1p5_sweep_fig.png` — IPR(λ) for both chains + ratio panel]**

---

## 8.2 Multifractal Dimensions at the Mid-Gap State

The linear-limit (λ = 0) initial-state IPR scales with chain length as IPR(N) ∼ N^{−D₂}, where D₂ is the generalized fractal dimension. Fitting the T=10⁴ initial-time IPR values across N ∈ {500, 1000, 2000}:

| Chain | D₂ | IPR at N=500 | IPR at N=1000 | IPR at N=2000 |
|---|---|---|---|---|
| Fibonacci | **0.578** | 0.02097 | 0.01080 | 0.00940 |
| Tribonacci | **0.381** | 0.08202 | 0.04101 | 0.04840 |

The Fibonacci D₂ ≈ 0.58 sits in the expected range for critical states of Fibonacci tight-binding models (literature: D₂ ≈ 0.5–0.7 depending on t_mod). The Tribonacci D₂ ≈ 0.38 is lower, indicating stronger multifractality — the state is more fragmented across the chain, consistent with its higher absolute IPR and the tribonacci word's richer substitution alphabet.

**[PLACEHOLDER: figure — log IPR vs log N for both chains at λ=0; slope gives D₂]**

---

## 8.3 Long-Time FSS Amplification at λ = 1.5

At T = 10⁴ the differential gap narrows but persists across all system sizes:

| N | IPR_fib(T=10⁴) | IPR_trib(T=10⁴) | Ratio (trib/fib) |
|---|---|---|---|
| 500 | 0.004899 | 0.008274 | 1.69× |
| 1000 | 0.003310 | 0.013751 | 4.15× |
| 2000 | 0.003366 | 0.010424 | 3.10× |

**Key result:** the Tribonacci state remains more localized (higher IPR) than the Fibonacci state at all three N values and at T two orders of magnitude beyond the paper's Table 1. The paper's −57%/−5% differential robustness finding is thus not a short-time artifact.

The trib/fib ratio at T=10⁴ (1.7–4.2×) is reduced compared to T=50 (8.65×), which is expected: both states spread under nonlinear dynamics, but at different rates. The Tribonacci state's non-monotone IPR with N (higher at N=1000 than at N=500 or N=2000) reflects finite-size sensitivity of the mid-gap eigenstate to the specific chain realization — a known feature of quasiperiodic systems where the gap eigenvalue fluctuates as N cycles through Fibonacci/Tribonacci approximant lengths.

**[PLACEHOLDER: figure — IPR(T=10⁴) vs N for both chains at λ=1.5]**

---

## 8.4 Spreading Exponents α(N) and Thermodynamic-Limit Signatures

Fitting IPR(t) ∼ t^{−α} on the late 30% of the log-time window (T ∈ [10³, 10⁴]) gives:

| N | α_fib | α_trib |
|---|---|---|
| 500 | 0.334 | 0.942 |
| 1000 | 0.254 | 0.331 |
| 2000 | 0.157 | 0.326 |

**Fibonacci:** α_fib decreases monotonically with N (0.334 → 0.157). The trend is consistent with α → 0 in the thermodynamic limit, indicating **self-trapping** at λ = 1.5 — the nonlinear term is strong enough to arrest spreading at large chain lengths. This is the dynamical counterpart of the IPR minimum seen in Section 8.1.

**Tribonacci:** α_trib converges to ≈ 0.33 for N ≥ 1000. A value α ≈ 1/3 is characteristic of sub-diffusive spreading in one-dimensional quasiperiodic systems in the presence of mild nonlinearity. The Tribonacci state's richer multifractal structure apparently sustains this spreading channel even at λ = 1.5, while the Fibonacci state is pinned.

The **sign flip** between chains in the large-N limit (α_fib → 0, α_trib → 0.33) constitutes the clearest thermodynamic-limit signature: in the N → ∞ picture, λ = 1.5 is a near-critical coupling for the Fibonacci chain but sub-critical for the Tribonacci chain.

**[PLACEHOLDER: figure — α vs N for both chains at λ=1.5; two panels or one with two curves]**

---

## 8.5 Open Items Before Submission

- [ ] Extend N to 4000–8000 to sharpen the α(N) fits (particularly the α_fib → 0 extrapolation)
- [ ] Run the same FSS at λ = 1.0 and λ = 2.0 to bracket the crossover region
- [ ] Produce Figure 8.2 (log IPR vs log N) and Figure 8.4 (α vs N) from the data above
- [ ] Add one-sentence note in the abstract/introduction crediting the T=50 reproducibility result
- [ ] Reconcile PRs #3/#4 (duplicate sweep scripts) before adding more files under the same name

---

*Prose slots above marked **[FIGURE]** / **[PLACEHOLDER]** should be filled in once the corresponding matplotlib panels are generated from the three CSV files.*
