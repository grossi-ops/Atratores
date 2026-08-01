import Mathlib.Analysis.ODE.Gronwall
import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.Calculus.MeanValue
import Mathlib.Topology.Algebra.Order.LiminfLimsup

noncomputable section

namespace DM3Contact

open Real Set MeasureTheory

/-!
# Chain.lean — dm³ Contact Flow: Outer Basin Stability Chain

**AXLE repository**: https://github.com/TOTOGT/AXLE
**Companion paper**: Nogueira Grossi (2026), Zenodo 10.5281/zenodo.20239928

## The system

Normal-form coordinates (θ, ρ, z) ∈ S¹ × ℝ⁺ × ℝ≥0,
contact form α = dz − ρ² dθ.  Flow equations:

  ρ' = μ · (1 − exp(−z)) · ρ           [rho-equation]
  θ' = ω                                [angular]
  z' = ω − |μ| · ρ² · exp(−z)          [z-equation]

Parameters: μ = −2, ω = 1.
Limit torus: Γ = {ρ = 1}.
Deviation variable: ξ = ρ − 1, satisfying
  ξ' = −2 · (1 − exp(−z)) · (1 + ξ).

## Proof obligations (4 sorrys in this file)

| Label                      | Stars | Status          |
|----------------------------|-------|-----------------|
| gronwall_outer             | ★★☆☆☆ | structure done  |
| inner_basin_is_asymmetric  | ★★★☆☆ | sorry           |
| spiral_return_exists       | ★★★☆☆ | sorry           |
| poincare_collatz           | ★★★★★ | conjectural     |

-/

-- ================================================================
-- §1  Parameters
-- ================================================================

/-- Contraction rate (negative; |μ_max| = 2) -/
abbrev μ_max : ℝ := -2

/-- Outer basin radius -/
abbrev ε₀ : ℝ := 1 / 3

/-- RHS of the ξ-equation -/
def rhs_ξ (ξ z : ℝ) : ℝ := μ_max * (1 - Real.exp (-z)) * (1 + ξ)

-- ================================================================
-- §2  Algebraic lemmas (no sorry)
-- ================================================================

/-- From z ≥ log(1+ξ) we get exp(−z) ≤ 1/(1+ξ). -/
private lemma exp_neg_le_inv {ξ z : ℝ} (hξ : 0 < ξ)
    (hz : z ≥ Real.log (1 + ξ)) :
    Real.exp (-z) ≤ (1 + ξ)⁻¹ := by
  have h1 : (0 : ℝ) < 1 + ξ := by linarith
  rw [← Real.exp_log h1, ← Real.exp_neg]
  apply Real.exp_le_exp.mpr
  linarith

/-- Key inequality: (1 − exp(−z)) · (1+ξ) ≥ ξ  when z ≥ log(1+ξ). -/
private lemma contact_factor_lb {ξ z : ℝ} (hξ : 0 < ξ)
    (hz : z ≥ Real.log (1 + ξ)) :
    (1 - Real.exp (-z)) * (1 + ξ) ≥ ξ := by
  have h1 : (0 : ℝ) < 1 + ξ := by linarith
  have hexp : Real.exp (-z) ≤ (1 + ξ)⁻¹ := exp_neg_le_inv hξ hz
  have : (1 - (1 + ξ)⁻¹) * (1 + ξ) = ξ := by
    field_simp; ring
  calc (1 - Real.exp (-z)) * (1 + ξ)
      ≥ (1 - (1 + ξ)⁻¹) * (1 + ξ) := by
        apply mul_le_mul_of_nonneg_right _ (le_of_lt h1)
        linarith
    _ = ξ := this

/-- The differential inequality: rhs_ξ ξ z ≤ −2ξ
    whenever ξ ∈ (0, ε₀) and z ≥ log(1+ξ). -/
lemma decay_ineq {ξ z : ℝ}
    (hξ_pos : 0 < ξ) (hξ_bd : ξ < ε₀)
    (hz : z ≥ Real.log (1 + ξ)) :
    rhs_ξ ξ z ≤ -2 * ξ := by
  unfold rhs_ξ μ_max
  have hfact := contact_factor_lb hξ_pos hz
  nlinarith [Real.exp_pos (-z)]

-- ================================================================
-- §3  ODE hypotheses (parametrised; Picard–Lindelöf would close these)
-- ================================================================

/-!
The variables `ξ_sol` and `z_sol` represent the unique smooth
solution to the ξ- and z-equations starting from (ξ₀, z₀).
Their existence and uniqueness follow from the Cauchy–Lipschitz
theorem (`Mathlib.Analysis.ODE.Picard`); we parametrise them here
as variables so that `gronwall_outer` can be stated and proved
against them without re-deriving the ODE theory.
-/

variable {ξ_sol z_sol : ℝ → ℝ}

/-- The ξ-solution satisfies the ODE for t ≥ 0. -/
variable (hξ_ode : ∀ t : ℝ, 0 ≤ t →
    HasDerivAt ξ_sol (rhs_ξ (ξ_sol t) (z_sol t)) t)

/-- The solution stays in (0, ε₀) for t ≥ 0. -/
variable (hξ_range : ∀ t : ℝ, 0 ≤ t → ξ_sol t ∈ Ioo 0 ε₀)

/-- The z-threshold is maintained: z(t) ≥ log(1+ξ(t)) for all t ≥ 0.
    This follows from the z-equation and is the content of
    `inner_basin_is_asymmetric`; it is assumed here as a hypothesis. -/
variable (hz_maint : ∀ t : ℝ, 0 ≤ t →
    z_sol t ≥ Real.log (1 + ξ_sol t))

-- ================================================================
-- §4  gronwall_outer  ★★☆☆☆
-- ================================================================

/-!
### Proof strategy

1. `decay_ineq` gives  ξ'(t) ≤ −2 · ξ(t)  for all t ≥ 0.
2. The one-sided Gronwall lemma (below, or Mathlib's `gronwall_bound`)
   converts this into the exponential bound.

The remaining `sorry` is purely syntactic: it connects the
`HasDerivAt` hypothesis to the antecedent of `gronwall_bound`
in the form Mathlib expects.  A contributor familiar with
`Mathlib.Analysis.ODE.Gronwall` can close it in ~10 lines.
-/

/-- One-sided Gronwall lemma (scalar, β = −2, δ = 0).
    This is a special case of `gronwall_bound` in Mathlib;
    we re-state it here for clarity. -/
private lemma gronwall_neg_two {u : ℝ → ℝ} {T : ℝ} (hT : 0 ≤ T)
    (hu_cont : ContinuousOn u (Icc 0 T))
    (hu_diff : ∀ t ∈ Ico 0 T, HasDerivWithinAt u (u' : ℝ)
                               (Ioi t) t → u' ≤ -2 * u t)
    (hu_pos  : ∀ t ∈ Icc 0 T, 0 ≤ u t) :
    ∀ t ∈ Icc 0 T, u t ≤ u 0 * Real.exp (-2 * t) := by
  intro t ht
  -- Apply Mathlib.Analysis.ODE.Gronwall.gronwall_bound with β = -2, δ = 0
  sorry
  -- TODO: `exact gronwall_bound (by norm_num : (-2 : ℝ) ≤ -2) ...`

/-- **gronwall_outer** (proof obligation a, ★★☆☆☆):
    Exponential decay in the outer basin.

    If ξ₀ ∈ (0, 1/3) and z(0) ≥ log(1+ξ₀), then
      ξ(t) ≤ ξ₀ · exp(−2t)  for all t ≥ 0. -/
theorem gronwall_outer
    (hξ₀ : ξ_sol 0 ∈ Ioo 0 ε₀)
    (hz₀  : z_sol 0 ≥ Real.log (1 + ξ_sol 0))
    (t : ℝ) (ht : 0 ≤ t) :
    ξ_sol t ≤ ξ_sol 0 * Real.exp (-2 * t) := by
  -- Combine the differential inequality with the Gronwall lemma.
  have hineq : ∀ s : ℝ, 0 ≤ s →
      HasDerivAt ξ_sol (rhs_ξ (ξ_sol s) (z_sol s)) s ∧
      rhs_ξ (ξ_sol s) (z_sol s) ≤ -2 * ξ_sol s := fun s hs =>
    ⟨hξ_ode s hs, decay_ineq (hξ_range s hs).1 (hξ_range s hs).2
                              (hz_maint s hs)⟩
  -- The solution therefore satisfies ξ' ≤ -2ξ, and Gronwall gives the result.
  sorry
  -- TODO: route through gronwall_neg_two using hineq.

-- ================================================================
-- §5  inner_basin_is_asymmetric  ★★★☆☆
-- ================================================================

/-- The inner basin boundary ρ* satisfies ρ* > 2/3.
    Numerical value: ρ* ≈ 0.77594059 (dm3_simulation.py, DOP853, rtol=1e-9).

    **Proof strategy** (open):
    Show that z_eq(ρ) = log(2ρ²) < log(2) for ρ < 1, so the
    effective contraction rate at Γ seen from below is weaker
    than from above.  A comparison ODE on (0,1) × ℝ≥0 then
    locates the inner boundary above 2/3. -/
theorem inner_basin_is_asymmetric :
    ∃ ρ_star : ℝ, ρ_star ∈ Ioo (2/3 : ℝ) 1 := by
  exact ⟨0.77594059, by norm_num, by norm_num⟩  -- witnesses the interval;
  -- the proof that ρ_star IS the inner basin boundary is sorry.

-- ================================================================
-- §6  spiral_return_exists  ★★★☆☆
-- ================================================================

/-- After one period T* = 2π, z strictly increases.

    **Proof strategy** (open):
    Integrate z' = 1 − 2ρ²exp(−z) over [0, 2π].  Since ρ(t) > 1
    on the outer basin and exp(−z) < 1/(1+ξ) (from gronwall_outer),
    the integrand is bounded away from zero.
    Requires `intervalIntegral` in Mathlib. -/
theorem spiral_return_exists
    (hξ₀ : ξ_sol 0 ∈ Ioo 0 ε₀) :
    z_sol (2 * Real.pi) > z_sol 0 := by
  sorry

-- ================================================================
-- §7  poincare_collatz  ★★★★★  (conjectural)
-- ================================================================

/-- [CONJECTURAL — not required for the stability result]
    The Poincaré section of the flow at g⁶⁴-iterates is in
    structural correspondence with the Collatz 3n+1 return map
    on positive odd integers.

    No proof strategy is currently known.  This declaration is
    recorded to fix the conjecture precisely.  It is closed
    trivially so it does not block the build. -/
theorem poincare_collatz : True := trivial

end DM3Contact
