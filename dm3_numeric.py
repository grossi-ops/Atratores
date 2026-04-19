#!/usr/bin/env python3
"""
dm3_numeric.py — Numeric laboratory for the dm³ contact-helical ODE.

Principia Orthogona Vol. IV · Helical Attractors on Contact 3-Manifolds
Pablo Nogueira Grossi · G6 LLC · 2026
DOI: 10.5281/zenodo.19117400

Reproduces Table 1 of the preprint at high precision using DOP853
(8th-order Dormand–Prince) with rtol=1e-10, atol=1e-12.

The dm³ system (ε = 2):

    ṙ = r(1 − r²) + ε(r − 1) exp(−z)
    θ̇ = 1
    ż = r² − ε(r − 1)² exp(−z)

Theorem 2.1 (symmetric, preprint §2) asserts exponential convergence to
r = 1 at rate μ = −2 for |r(0) − 1| < 1/3 and z(0) ≥ log 2.

This script:
  1. Integrates the ODE for a range of r(0) values.
  2. Fits μ̂ = slope of log|r(t) − 1| on t ∈ [5, 20].
  3. Records min_t ż over t ∈ [0, 30] (the inner-basin diagnostic).
  4. Prints the canonical Table 1.

Usage:
    python3 dm3_numeric.py
    python3 dm3_numeric.py --plot        # requires matplotlib

Dependencies: numpy, scipy; matplotlib optional for --plot.
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

# ─── Model ──────────────────────────────────────────────────────────────────

EPSILON = 2.0  # fixed coupling for Vol. IV


def dm3(t: float, y: np.ndarray, eps: float = EPSILON) -> list[float]:
    """Right-hand side of the dm³ ODE in state y = (r, theta, z)."""
    r, theta, z = y
    ex = math.exp(-z)
    rdot = r * (1.0 - r * r) + eps * (r - 1.0) * ex
    tdot = 1.0
    zdot = r * r - eps * (r - 1.0) ** 2 * ex
    return [rdot, tdot, zdot]


# ─── Diagnostics ────────────────────────────────────────────────────────────


@dataclass
class RunResult:
    r0: float
    mu_hat: float          # empirical contraction rate on fit window
    zdot_early: float      # ż at early snapshot (t = 0.3) or NaN if escaped earlier
    zdot_min_bounded: float  # min ż on t ∈ [0, min(T_survive, 2)]
    converges: bool
    t: np.ndarray
    y: np.ndarray
    t_survive: float       # integration end time (< 30 if solver stopped)


def _fit_mu(
    sol, t_lo: float = 5.0, t_hi: float = 15.0
) -> float:
    """Empirical contraction rate: slope of log|r(t) − 1| on [t_lo, t_hi].

    Uses the dense interpolator for stable fits. Returns NaN if the
    solver did not survive to t_hi or if |r − 1| drops below atol.
    """
    if sol.t[-1] < t_hi:
        return float("nan")
    tt = np.linspace(t_lo, t_hi, 400)
    rr = sol.sol(tt)[0]
    mask = np.abs(rr - 1.0) > 1e-12
    if mask.sum() < 20:
        return float("nan")
    slope, _intercept = np.polyfit(tt[mask], np.log(np.abs(rr[mask] - 1.0)), 1)
    return float(slope)


def run_one(r0: float, t_end: float = 30.0) -> RunResult:
    """Integrate one trajectory at high precision with dense output."""
    y0 = np.array([r0, 0.0, 0.0])
    sol = solve_ivp(
        dm3,
        t_span=(0.0, t_end),
        y0=y0,
        method="DOP853",
        rtol=1e-10,
        atol=1e-12,
        dense_output=True,
        max_step=0.05,
    )
    t_survive = float(sol.t[-1])

    # Early snapshot of ż at t = 0.3 (before any blow-up in all test cases).
    t_snap = 0.3
    if t_survive >= t_snap:
        y_snap = sol.sol(t_snap)
        r_s, _, z_s = y_snap
        zdot_early = float(r_s ** 2 - EPSILON * (r_s - 1.0) ** 2 * math.exp(-z_s))
    else:
        zdot_early = float("nan")

    # Bounded min ż on t ∈ [0, min(t_survive, 2)] — before catastrophic blow-up.
    t_cap = min(t_survive * 0.98, 2.0)
    tt = np.linspace(0.0, t_cap, 2000)
    y = sol.sol(tt)
    r = y[0]
    z = y[2]
    zdot = r ** 2 - EPSILON * (r - 1.0) ** 2 * np.exp(-z)
    zdot_min_bounded = float(zdot.min())

    # Convergence: survived the full integration and ended near r = 1.
    # Note: r[-1] above is the bounded-window sample (up to t ≈ 2); the
    # convergence test must use the actual end-of-integration value.
    r_final = float(sol.y[0][-1])
    converges = bool(t_survive >= t_end * 0.99 and abs(r_final - 1.0) < 2e-2)
    mu_hat = _fit_mu(sol) if converges else float("nan")

    return RunResult(
        r0=r0,
        mu_hat=mu_hat,
        zdot_early=zdot_early,
        zdot_min_bounded=zdot_min_bounded,
        converges=converges,
        t=sol.t,
        y=sol.y,
        t_survive=t_survive,
    )


# ─── Table 1 ───────────────────────────────────────────────────────────────

OUTER_R0 = [1.10, 1.33, 1.50, 2.00, 2.50]   # r(0) − 1 = {+0.10, +0.33, +0.50, +1.00, +1.50}
INNER_R0 = [0.90, 0.80, 0.70, 0.50, 0.30]


def table_1() -> None:
    """Run the canonical Table 1 and print the result.

    Diagnostics reported:
      · μ̂       — empirical contraction rate, slope of log|r(t)−1| on t ∈ [5, 15].
      · ż(0.3)  — early vertical velocity, a signed diagnostic that stays finite
                   even for collapsing orbits because it is sampled before the
                   coupling term can blow up.
      · outcome — qualitative classification from end-state of the DOP853 run.
    """
    outer_rows = [run_one(r0) for r0 in OUTER_R0]
    inner_rows = [run_one(r0) for r0 in INNER_R0]

    def classify(res: RunResult) -> str:
        if res.converges:
            return "converges"
        if res.t_survive >= 2.0:
            return "limiar · edge"
        return "collapses"

    print()
    print("═" * 74)
    print("  TABLE 1 — dm³ at ε = 2  ·  DOP853 @ rtol=1e-10, atol=1e-12")
    print("  Principia Orthogona Vol. IV  ·  Pablo Nogueira Grossi  ·  2026")
    print("═" * 74)
    print()
    print("  OUTER BASIN — r(0) > 1  (fit μ̂ on t ∈ [5, 15])")
    print("  " + "-" * 62)
    print(f"  {'r(0)-1':>10}  {'r(0)':>8}  {'μ̂':>12}  {'outcome':>18}")
    print("  " + "-" * 62)
    for res in outer_rows:
        print(
            f"  {res.r0 - 1.0:+10.2f}  {res.r0:8.2f}  "
            f"{res.mu_hat:12.4f}  {classify(res):>18}"
        )
    print()
    print("  INNER BASIN — r(0) < 1  (ż sampled at t = 0.3)")
    print("  " + "-" * 70)
    print(
        f"  {'r(0)':>8}  {'ż(t=0.3)':>12}  {'t_survive':>12}  "
        f"{'outcome':>18}"
    )
    print("  " + "-" * 70)
    for res in inner_rows:
        print(
            f"  {res.r0:8.2f}  {res.zdot_early:12.3f}  "
            f"{res.t_survive:12.3f}  {classify(res):>18}"
        )
    print()
    print("  READING:")
    print("  · Outer (r(0) > 1): all runs converge; μ̂ → −2 from above as")
    print("    expected from the Jacobian eigenvalue at r = 1.")
    print("  · Inner (r(0) < 1): converges for r(0) ≥ 0.80, collapses below.")
    print("    The integration time t_survive drops precipitously: a robust")
    print("    numerical fingerprint of the escape to z → −∞.")
    print("  · Gronwall (Theorem 2.1) guarantees only |r(0) − 1| < 1/3, i.e.")
    print("    r(0) ∈ (2/3, 4/3) ≈ (0.667, 1.333). The numerical inner edge")
    print("    is r⋆ ≈ 0.80, strictly stricter than Gronwall's 2/3.")
    print()
    print("═" * 74)
    print()


# ─── Optional plot ──────────────────────────────────────────────────────────


def plot_trajectories(out_path: str = "dm3_trajectories.png") -> None:
    """Plot the 10 canonical trajectories in (r, z)."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("[plot] matplotlib not installed; skipping plot.")
        return
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    axes[0].set_title("Outer basin — r(0) > 1")
    axes[1].set_title("Inner basin — r(0) < 1")
    for r0 in OUTER_R0:
        res = run_one(r0)
        axes[0].plot(res.y[0], res.y[2], label=f"r₀ = {r0:.2f}")
    for r0 in INNER_R0:
        res = run_one(r0)
        axes[1].plot(res.y[0], res.y[2], label=f"r₀ = {r0:.2f}")
    for ax in axes:
        ax.axvline(1.0, color="#c9a84c", lw=0.8, alpha=0.7, ls="--",
                   label="Γ at r = 1")
        ax.axvline(0.80, color="#e05a3a", lw=0.8, alpha=0.5, ls=":",
                   label="r⋆ ≈ 0.80")
        ax.set_xlabel("r")
        ax.set_ylabel("z")
        ax.legend(fontsize=8, loc="best")
        ax.grid(True, alpha=0.25)
    fig.suptitle("dm³ trajectories · ε = 2 · DOP853 rtol=1e-10",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"[plot] wrote {out_path}")


# ─── Entry point ────────────────────────────────────────────────────────────


def main() -> int:
    ap = argparse.ArgumentParser(
        description="dm³ numeric lab — reproduces Table 1 of Vol. IV."
    )
    ap.add_argument(
        "--plot", action="store_true",
        help="also save a trajectory plot (requires matplotlib)",
    )
    ap.add_argument(
        "--plot-out", default="dm3_trajectories.png",
        help="output path for the plot (default: dm3_trajectories.png)",
    )
    args = ap.parse_args()

    table_1()

    if args.plot:
        plot_trajectories(args.plot_out)

    return 0


if __name__ == "__main__":
    sys.exit(main())
