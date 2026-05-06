#!/usr/bin/env python3
"""
analyze_long_time.py
====================
Analysis companion for `dnls_long_time.py`.

Reads `ipr_vs_time.csv` (long-format: time, lambda, chain, IPR, norm) and
produces:

  1. A sanity report (text) covering:
       - lambda = 0 flatness check  (linear-limit baseline)
       - norm-conservation summary  (DOP853 drift per run)
       - t = 0 IPR values           (tie-back to Table 1 of the paper)
  2. fig_long_ipr_vs_t.png         - IPR(t) on log-log, both chains, all lambdas
  3. fig_long_lambda0_check.png    - flatness of the linear-limit run
  4. fig_long_alpha_fit.png        - late-time fits of IPR(t) ~ t^(-alpha)
  5. spreading_exponents.csv       - fitted alpha per (chain, lambda)

``--report`` mode (Brief 4 / T=10^6 verification)
--------------------------------------------------
  python3 analyze_long_time.py --report --csv data/ipr_lambda1p5_N1000_T1e6.csv \
      --csv-t5 data/ipr_lambda1p5_N1000_T1e5.csv \
      --fig figures/T1e6_lambda1p5_N1000.png

  Produces verbatim sections [1]-[5] as printed output, plus the figure.
  fit_alpha_full() is used internally for three time windows:
      t > 1e4, t > 1e5, t > 3e5.

Companion to:
  "Differential Nonlinear Robustness of Critical States in Fibonacci and
   Tribonacci Substitution Chains"
  Pablo Nogueira Grossi, G6 LLC (2026)
  DOI: 10.5281/zenodo.20026943

Usage
-----
    python3 analyze_long_time.py                          # uses ipr_vs_time.csv
    python3 analyze_long_time.py --csv long_run.csv       # custom input
    python3 analyze_long_time.py --report \
        --csv data/ipr_lambda1p5_N1000_T1e6.csv \
        --csv-t5 data/ipr_lambda1p5_N1000_T1e5.csv \
        --fig figures/T1e6_lambda1p5_N1000.png

Author
------
    Pablo Nogueira Grossi  |  ORCID: 0009-0000-6496-2186
    G6 LLC, Newark NJ  |  GitHub: https://github.com/TOTOGT/AXLE

License: MIT
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# 1. CSV loader
# ---------------------------------------------------------------------------

def load_csv(path: str) -> dict[tuple[str, float], dict[str, np.ndarray]]:
    """
    Load the long-format CSV into a {(chain, lambda): {time, IPR, norm}} dict.
    """
    rows_by_key: dict[tuple[str, float], list[tuple[float, float, float]]] = (
        defaultdict(list)
    )
    with open(path, "r", newline="") as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            key = (r["chain"], float(r["lambda"]))
            rows_by_key[key].append(
                (float(r["time"]), float(r["IPR"]), float(r["norm"]))
            )

    out: dict[tuple[str, float], dict[str, np.ndarray]] = {}
    for key, rows in rows_by_key.items():
        rows.sort(key=lambda x: x[0])
        t = np.array([r[0] for r in rows])
        ipr = np.array([r[1] for r in rows])
        norm = np.array([r[2] for r in rows])
        out[key] = {"time": t, "IPR": ipr, "norm": norm}
    return out


# ---------------------------------------------------------------------------
# 2. Sanity report
# ---------------------------------------------------------------------------

def sanity_report(data: dict[tuple[str, float], dict[str, np.ndarray]]) -> dict:
    """
    Print and return a structured sanity report.

    Hard checks:
      - lambda=0 IPR must be flat to ~5 sig figs (linear-limit eigenstate)
      - norm drift max should be << 1e-3 for any reasonable run
    """
    lambdas = sorted({lam for _, lam in data.keys()})
    chains = sorted({ch for ch, _ in data.keys()})

    print("=" * 72)
    print("SANITY REPORT")
    print("=" * 72)

    # --- lambda=0 flatness check ---
    print("\n[1] Linear-limit (lambda=0) flatness check")
    print("    IPR(t) at lambda=0 should be constant -- the t=0 eigenstate")
    print("    evolves only in phase under linear dynamics.")
    print()
    print(f"    {'chain':>12} {'IPR(t=0)':>14} {'IPR(t=T)':>14} "
          f"{'max|dIPR|':>14} {'rel.var.':>12}")
    print("    " + "-" * 68)
    flatness = {}
    for ch in chains:
        key = (ch, 0.0)
        if key not in data:
            print(f"    {ch:>12}  no lambda=0 row -- linear-limit check skipped")
            continue
        ipr = data[key]["IPR"]
        ipr0 = float(ipr[0])
        iprT = float(ipr[-1])
        max_dev = float(np.max(np.abs(ipr - ipr0)))
        rel = max_dev / abs(ipr0) if ipr0 != 0 else float("inf")
        flatness[ch] = (ipr0, iprT, max_dev, rel)
        print(f"    {ch:>12} {ipr0:>14.8f} {iprT:>14.8f} "
              f"{max_dev:>14.2e} {rel:>12.2e}")

    # --- t=0 IPR (tie-back to Table 1) ---
    print("\n[2] t = 0 IPR values  (tie-back to Table 1 of the paper)")
    print("    These should match the lambda=0 column of Table 1 exactly.")
    print()
    print(f"    {'chain':>12} {'IPR(0)':>14}")
    print("    " + "-" * 28)
    for ch in chains:
        # any lambda will do; t=0 is identical (mid-gap eigenstate)
        for lam in lambdas:
            key = (ch, lam)
            if key in data:
                print(f"    {ch:>12} {float(data[key]['IPR'][0]):>14.8f}")
                break

    # --- norm-conservation summary ---
    print("\n[3] Norm conservation across all runs")
    print("    DOP853 at rtol=1e-8 should keep |norm-1| << 1e-5 over T=10^3.")
    print()
    print(f"    {'chain':>12} {'lambda':>8} {'max|norm-1|':>14}")
    print("    " + "-" * 36)
    norm_summary = []
    worst_drift = 0.0
    for ch in chains:
        for lam in lambdas:
            key = (ch, lam)
            if key not in data:
                continue
            drift = float(np.max(np.abs(data[key]["norm"] - 1.0)))
            worst_drift = max(worst_drift, drift)
            tag = "" if drift < 1e-5 else "  <-- check"
            print(f"    {ch:>12} {lam:>8.2f} {drift:>14.2e}{tag}")
            norm_summary.append((ch, lam, drift))

    # --- final IPR vs lambda summary (compact table) ---
    print("\n[4] Final-time IPR(T) by chain and lambda  (the long-time headline)")
    print()
    header = "    " + f"{'lambda':>8}"
    for ch in chains:
        header += f" {('IPR_' + ch[:4]):>14}"
    if "fibonacci" in chains and "tribonacci" in chains:
        header += f" {'ratio_t/f':>12}"
    print(header)
    print("    " + "-" * (len(header) - 4))

    for lam in lambdas:
        line = f"    {lam:>8.2f}"
        ipr_by_chain = {}
        for ch in chains:
            key = (ch, lam)
            if key in data:
                v = float(data[key]["IPR"][-1])
                ipr_by_chain[ch] = v
                line += f" {v:>14.8f}"
            else:
                line += f" {'-':>14}"
        if "fibonacci" in ipr_by_chain and "tribonacci" in ipr_by_chain:
            line += f" {ipr_by_chain['tribonacci']/ipr_by_chain['fibonacci']:>12.4f}"
        print(line)

    print()
    print("=" * 72)
    return {
        "flatness": flatness,
        "norm_summary": norm_summary,
        "worst_drift": worst_drift,
    }


# ---------------------------------------------------------------------------
# 3. Plots
# ---------------------------------------------------------------------------

def plot_ipr_vs_t(
    data: dict[tuple[str, float], dict[str, np.ndarray]],
    out_path: str = "fig_long_ipr_vs_t.png",
) -> None:
    """
    Log-log IPR(t) for both chains, all lambdas. Solid = tribonacci, dashed = fibonacci.
    """
    chains = sorted({ch for ch, _ in data.keys()})
    lambdas = sorted({lam for _, lam in data.keys()})

    fig, ax = plt.subplots(figsize=(8, 6), dpi=140)
    cmap = plt.cm.viridis(np.linspace(0.05, 0.95, len(lambdas)))
    style = {"fibonacci": "--", "tribonacci": "-"}
    width = {"fibonacci": 1.3, "tribonacci": 1.8}

    for ch in chains:
        for i, lam in enumerate(lambdas):
            key = (ch, lam)
            if key not in data:
                continue
            t = data[key]["time"]
            ipr = data[key]["IPR"]
            mask = t > 0
            ax.loglog(
                t[mask],
                ipr[mask],
                style.get(ch, "-"),
                color=cmap[i],
                lw=width.get(ch, 1.5),
                label=f"{ch[:4]}  lambda={lam:.1f}",
            )

    ax.set_xlabel("time  t")
    ax.set_ylabel("IPR(t)")
    ax.set_title("Long-time IPR on Fibonacci (dashed) and Tribonacci (solid) chains")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(ncol=2, fontsize=8, loc="best")
    plt.tight_layout()
    plt.savefig(out_path, dpi=140)
    plt.close()
    print(f"  -> {out_path}")


def plot_lambda0_check(
    data: dict[tuple[str, float], dict[str, np.ndarray]],
    out_path: str = "fig_long_lambda0_check.png",
) -> None:
    """
    Linear-limit sanity plot: IPR(t) at lambda=0 should be visually flat.
    """
    chains = sorted({ch for ch, _ in data.keys()})
    fig, ax = plt.subplots(figsize=(8, 5), dpi=140)
    for ch in chains:
        key = (ch, 0.0)
        if key not in data:
            continue
        t = data[key]["time"]
        ipr = data[key]["IPR"]
        ax.semilogx(t[t > 0], ipr[t > 0], lw=1.6, label=f"{ch}  IPR(0)={ipr[0]:.5f}")
    ax.set_xlabel("time  t")
    ax.set_ylabel("IPR(t)  at lambda = 0")
    ax.set_title("Linear-limit sanity check (IPR should be flat)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    plt.tight_layout()
    plt.savefig(out_path, dpi=140)
    plt.close()
    print(f"  -> {out_path}")


# ---------------------------------------------------------------------------
# 4. Spreading-exponent fit  (open question 3)
# ---------------------------------------------------------------------------

def fit_alpha(t: np.ndarray, ipr: np.ndarray, late_frac: float = 0.3
              ) -> tuple[float, float, int]:
    """
    Fit IPR(t) ~ A * t^(-alpha) on the late-time tail.

    Strategy: take the last `late_frac` of points by log(t), do a linear
    regression of log(IPR) on log(t). Returns (alpha, log_A, n_points).

    A positive alpha means the participation length L_eff ~ 1 / IPR is
    spreading; alpha = 1 is ballistic / fully extended; alpha < 1 is
    sub-diffusive; alpha approx 0 is self-trapping.
    """
    mask = (t > 0) & (ipr > 0)
    t = t[mask]
    ipr = ipr[mask]
    if len(t) < 6:
        return float("nan"), float("nan"), 0

    log_t = np.log(t)
    log_ipr = np.log(ipr)
    cutoff = log_t[0] + (1.0 - late_frac) * (log_t[-1] - log_t[0])
    sel = log_t >= cutoff
    if sel.sum() < 4:
        return float("nan"), float("nan"), int(sel.sum())

    slope, intercept = np.polyfit(log_t[sel], log_ipr[sel], 1)
    return float(-slope), float(intercept), int(sel.sum())


def fit_alpha_full(
    t: np.ndarray,
    ipr: np.ndarray,
    t_min: float,
) -> tuple[float, float, float, float, int]:
    """
    Fit IPR(t) ~ A * t^(-alpha) for all points with t > t_min.

    Uses OLS on log(IPR) vs log(t) with the standard stderr on the slope.

    Parameters
    ----------
    t, ipr : arrays of times and IPR values
    t_min  : lower time boundary for the fit window

    Returns
    -------
    alpha      : spreading exponent (positive = spreading)
    alpha_se   : standard error on alpha from OLS residuals
    log_A      : ln(amplitude A) from the fit
    r2         : coefficient of determination R^2
    n_pts      : number of points used
    """
    mask = (t > t_min) & (ipr > 0)
    n = int(mask.sum())
    if n < 4:
        return float("nan"), float("nan"), float("nan"), float("nan"), n

    lt = np.log(t[mask])
    li = np.log(ipr[mask])

    # OLS: li = intercept + slope * lt
    lt_bar = lt.mean()
    li_bar = li.mean()
    Sxx = float(np.sum((lt - lt_bar) ** 2))
    Sxy = float(np.sum((lt - lt_bar) * (li - li_bar)))
    if Sxx == 0.0:
        return float("nan"), float("nan"), float("nan"), float("nan"), n

    slope = Sxy / Sxx
    intercept = li_bar - slope * lt_bar

    # R^2
    li_pred = intercept + slope * lt
    ss_res = float(np.sum((li - li_pred) ** 2))
    ss_tot = float(np.sum((li - li_bar) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    # stderr on slope: se(slope) = sqrt(ss_res / ((n-2) * Sxx))
    if n > 2:
        se_slope = float(np.sqrt(ss_res / ((n - 2) * Sxx)))
    else:
        se_slope = float("nan")

    return float(-slope), se_slope, float(intercept), r2, n


# ---------------------------------------------------------------------------
# 5. Brief 4 verification report  (T=10^6, N=1000, λ=1.5)
# ---------------------------------------------------------------------------

_WINDOWS = [
    ("t > 1e4",  1e4,  "t > 10⁴"),
    ("t > 1e5",  1e5,  "t > 10⁵"),
    ("t > 3e5",  3e5,  "t > 3×10⁵"),
]


def _ipr_at_t(data_entry: dict, t_target: float) -> float:
    """Return IPR value at the checkpoint nearest to t_target."""
    t = data_entry["time"]
    ipr = data_entry["IPR"]
    idx = int(np.argmin(np.abs(t - t_target)))
    return float(ipr[idx])


def verification_report(
    data6: dict[tuple[str, float], dict[str, np.ndarray]],
    data5: dict[tuple[str, float], dict[str, np.ndarray]] | None = None,
    lam: float = 1.5,
) -> dict:
    """
    Print verbatim Brief 4 sections [1]–[4] and return a results dict.

    Parameters
    ----------
    data6 : loaded CSV for T=10^6 run
    data5 : loaded CSV for T=10^5 run (used for T=10^5 baseline column)
    lam   : lambda value to analyse (default 1.5)
    """
    chains = ["fibonacci", "tribonacci"]

    # ---- [1] Final-time IPR comparison -----------------------------------
    print("=" * 72)
    print(f"[1] Final-time IPR comparison  (λ={lam}, N=1000)")
    print("=" * 72)

    # Collect IPR at three time points from data6 (which covers 1→10^6)
    # T=10^4 and T=10^5 snapshots come from data6 checkpoints;
    # T=10^5 cross-check comes from data5 if available.
    ipr_t4 = {}
    ipr_t5 = {}
    ipr_t6 = {}
    norm_info = {}

    for ch in chains:
        key = (ch, lam)
        if key in data6:
            d = data6[key]
            ipr_t4[ch] = _ipr_at_t(d, 1e4)
            ipr_t6[ch] = float(d["IPR"][-1])
            norm_info[ch] = float(np.max(np.abs(d["norm"] - 1.0)))
        if data5 is not None and key in data5:
            ipr_t5[ch] = float(data5[key]["IPR"][-1])
        elif key in data6:
            ipr_t5[ch] = _ipr_at_t(data6[key], 1e5)

    def ratio_str(d: dict, a: str, b: str) -> str:
        if a in d and b in d and d[b] != 0:
            return f"{d[a]/d[b]:.4f}"
        return "n/a"

    print(f"\n  λ={lam}, N=1000:")
    for label, ipr_d in [("T=10⁴", ipr_t4), ("T=10⁵", ipr_t5), ("T=10⁶", ipr_t6)]:
        fib = ipr_d.get("fibonacci", float("nan"))
        trib = ipr_d.get("tribonacci", float("nan"))
        r_ft = fib / trib if trib != 0 else float("nan")
        r_tf = trib / fib if fib != 0 else float("nan")
        print(
            f"  {label} : IPR_fib = {fib:.6f}, IPR_trib = {trib:.6f}, "
            f"ratio fib/trib = {r_ft:.4f}, trib/fib = {r_tf:.4f}×"
        )

    # Cross-over check: scan data6 for first t where trib/fib < 1
    crossover_t = None
    if ("fibonacci", lam) in data6 and ("tribonacci", lam) in data6:
        t_f = data6[("fibonacci", lam)]["time"]
        ipr_f = data6[("fibonacci", lam)]["IPR"]
        ipr_t = data6[("tribonacci", lam)]["IPR"]
        # Interpolate onto the fib time grid using tribonacci checkpoints
        t_t = data6[("tribonacci", lam)]["time"]
        ipr_t_interp = np.interp(t_f, t_t, ipr_t)
        ratio_series = ipr_t_interp / ipr_f
        cross_mask = ratio_series < 1.0
        if cross_mask.any():
            crossover_t = float(t_f[cross_mask][0])
            print(f"\n  Cross-over (trib/fib < 1):  YES at t ≈ {crossover_t:.3e}")
        else:
            trib_fib_final = ipr_t6.get("tribonacci", float("nan")) / ipr_t6.get("fibonacci", float("nan"))
            print(f"\n  Cross-over (trib/fib < 1):  NO (trib/fib = {trib_fib_final:.4f} at T=10⁶)")

    # ---- [2] α-fit table at three windows --------------------------------
    print()
    print("=" * 72)
    print("[2] α-fit table at three time windows  (IPR ~ t^{-α})")
    print("=" * 72)
    print(f"\n  {'window':<12} {'chain':<12} {'α':>8} {'α_stderr':>10} {'R²':>8} {'n_pts':>6}")
    print("  " + "-" * 56)

    fits: dict[tuple[str, str], tuple[float, float, float, float, int]] = {}
    for win_key, t_min, win_label in _WINDOWS:
        for ch in chains:
            key = (ch, lam)
            if key not in data6:
                continue
            t = data6[key]["time"]
            ipr = data6[key]["IPR"]
            alpha, alpha_se, log_A, r2, n = fit_alpha_full(t, ipr, t_min)
            fits[(win_key, ch)] = (alpha, alpha_se, log_A, r2, n)
            a_str = f"{alpha:.4f}" if np.isfinite(alpha) else "  n/a "
            se_str = f"{alpha_se:.4f}" if np.isfinite(alpha_se) else "  n/a "
            r2_str = f"{r2:.4f}" if np.isfinite(r2) else "  n/a "
            print(f"  {win_label:<12} {ch:<12} {a_str:>8} {se_str:>10} {r2_str:>8} {n:>6d}")
        print()

    # ---- [3] Cross-over prediction from t > 10^5 window -----------------
    print("=" * 72)
    print("[3] Cross-over prediction  (from t > 10⁵ window)")
    print("=" * 72)

    fib6 = ipr_t6.get("fibonacci", float("nan"))
    trib6 = ipr_t6.get("tribonacci", float("nan"))
    ratio_now = trib6 / fib6 if fib6 != 0 else float("nan")

    win5_key = "t > 1e5"
    a_fib = fits.get((win5_key, "fibonacci"), (float("nan"),) * 5)[0]
    a_trib = fits.get((win5_key, "tribonacci"), (float("nan"),) * 5)[0]
    d_alpha = a_trib - a_fib

    if np.isfinite(d_alpha) and d_alpha > 0 and np.isfinite(ratio_now) and ratio_now > 1.0:
        t_cross = 1e6 * ratio_now ** (1.0 / d_alpha)
        brief1_consistent = abs(np.log10(t_cross) - np.log10(1.4e6)) < 0.5
        print(f"\n  trib/fib at T=10⁶        = {ratio_now:.4f}")
        print(f"  α_trib (t > 10⁵)         = {a_trib:.4f}")
        print(f"  α_fib  (t > 10⁵)         = {a_fib:.4f}")
        print(f"  Δα = α_trib − α_fib       = {d_alpha:.4f}")
        print(f"\n  t_cross_predicted = 10⁶ × ({ratio_now:.4f})^(1/{d_alpha:.4f})")
        print(f"                    = {t_cross:.3e}")
        if brief1_consistent:
            print(f"  → Consistent with Brief 1 prediction (~1.4×10⁶).")
        else:
            print(f"  → Significantly different from Brief 1 prediction (~1.4×10⁶).")
    elif np.isfinite(d_alpha) and d_alpha <= 0:
        print(f"\n  α_trib ({a_trib:.4f}) ≤ α_fib ({a_fib:.4f}): Δα ≤ 0.")
        print(f"  No power-law cross-over predicted under this extrapolation.")
    elif crossover_t is not None:
        print(f"\n  Cross-over already occurred within T=10⁶ at t ≈ {crossover_t:.3e}.")
    else:
        print(f"\n  trib/fib at T=10⁶ = {ratio_now:.4f}")
        print(f"  α data insufficient for cross-over prediction.")

    # ---- [4] Norm conservation -------------------------------------------
    print()
    print("=" * 72)
    print("[4] Norm conservation")
    print("=" * 72)
    print(f"\n  {'chain':<12} {'max|norm-1|':>14}")
    print("  " + "-" * 28)
    max_drift_all = 0.0
    for ch in chains:
        key = (ch, lam)
        if key not in data6:
            continue
        drift = norm_info.get(ch, float("nan"))
        if np.isfinite(drift):
            max_drift_all = max(max_drift_all, drift)
        flag = "  ← EXCEEDS 1e-4" if drift > 1e-4 else ""
        print(f"  {ch:<12} {drift:>14.2e}{flag}")
    print(f"\n  Max drift overall: {max_drift_all:.2e}")
    if max_drift_all > 1e-4:
        print("  WARNING: norm leak exceeds 1e-4 threshold.")
    else:
        print("  Norm conservation acceptable (< 1e-4).")

    print()
    print("=" * 72)

    return {
        "ipr_t4": ipr_t4,
        "ipr_t5": ipr_t5,
        "ipr_t6": ipr_t6,
        "fits": fits,
        "crossover_t": crossover_t,
        "norm_info": norm_info,
    }


def plot_verification(
    data6: dict[tuple[str, float], dict[str, np.ndarray]],
    fits: dict,
    crossover_t: float | None,
    ipr_t6: dict,
    lam: float = 1.5,
    fig_out: str = "figures/T1e6_lambda1p5_N1000.png",
) -> None:
    """
    Two-panel figure for Brief 4:
      Top   : log-log IPR(t) for both chains + α-fit lines for each window
      Bottom: ratio trib/fib vs t (log-x, linear-y) with reference lines
    """
    import os
    os.makedirs(os.path.dirname(fig_out) or ".", exist_ok=True)

    chains = ["fibonacci", "tribonacci"]
    colors = {"fibonacci": "#1f77b4", "tribonacci": "#d62728"}
    lstyles = {"fibonacci": "--", "tribonacci": "-"}

    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(9, 10), dpi=140)

    # -- top panel: IPR(t) --------------------------------------------------
    for ch in chains:
        key = (ch, lam)
        if key not in data6:
            continue
        t = data6[key]["time"]
        ipr = data6[key]["IPR"]
        mask = (t > 0) & (ipr > 0)
        ax_top.loglog(
            t[mask], ipr[mask],
            lstyles[ch], color=colors[ch], lw=2.0,
            label=f"{ch}",
        )
        # Draw α-fit lines for each window
        win_styles = [("t > 1e4", (5, 2), 1.4), ("t > 1e5", (3, 1, 1, 1), 1.9), ("t > 3e5", (1, 1), 2.2)]
        for win_key, dash, lw in win_styles:
            res = fits.get((win_key, ch))
            if res is None:
                continue
            alpha, _, log_A, r2, n = res
            if not np.isfinite(alpha) or n < 4:
                continue
            t_min = {"t > 1e4": 1e4, "t > 1e5": 1e5, "t > 3e5": 3e5}[win_key]
            t_fit = t[mask & (t >= t_min)]
            if len(t_fit) < 2:
                continue
            fit_line = np.exp(log_A) * t_fit ** (-alpha)
            ax_top.loglog(
                t_fit, fit_line,
                color=colors[ch], lw=lw, dashes=dash,
                label=f"  {ch[:4]} α={alpha:.3f} ({win_key})",
            )

    ax_top.set_xlabel("time  t")
    ax_top.set_ylabel("IPR(t)")
    ax_top.set_title(f"T=10⁶ verification: IPR(t) at λ={lam}, N=1000\n"
                     "(solid/dashed = chain; dotted = α-fit lines)")
    ax_top.grid(True, which="both", alpha=0.3)
    ax_top.legend(fontsize=7, ncol=2, loc="lower left")

    # -- bottom panel: ratio trib/fib vs t ----------------------------------
    key_f = ("fibonacci", lam)
    key_t = ("tribonacci", lam)
    if key_f in data6 and key_t in data6:
        t_f = data6[key_f]["time"]
        ipr_f = data6[key_f]["IPR"]
        t_t = data6[key_t]["time"]
        ipr_t_vals = data6[key_t]["IPR"]
        # Interpolate tribonacci onto fibonacci time grid
        ipr_t_interp = np.interp(t_f, t_t, ipr_t_vals)
        ratio = ipr_t_interp / np.where(ipr_f > 0, ipr_f, np.nan)
        mask_r = (t_f > 0) & np.isfinite(ratio) & (ratio > 0)
        ax_bot.semilogx(t_f[mask_r], ratio[mask_r], "k-", lw=2.0, label="trib/fib")
        ax_bot.axhline(1.0, color="gray", lw=1.5, ls="--", label="ratio = 1 (cross-over)")
        if crossover_t is not None:
            ax_bot.axvline(crossover_t, color="crimson", lw=1.5, ls=":", label=f"cross-over t={crossover_t:.1e}")
        else:
            # draw predicted cross-over from t>1e5 window
            a_fib5 = fits.get(("t > 1e5", "fibonacci"), (float("nan"),) * 5)[0]
            a_trib5 = fits.get(("t > 1e5", "tribonacci"), (float("nan"),) * 5)[0]
            ratio_now = ipr_t6.get("tribonacci", float("nan")) / ipr_t6.get("fibonacci", float("nan"))
            d_alpha = a_trib5 - a_fib5
            if np.isfinite(d_alpha) and d_alpha > 0 and np.isfinite(ratio_now) and ratio_now > 1.0:
                t_pred = 1e6 * ratio_now ** (1.0 / d_alpha)
                if t_pred < 1e10:
                    ax_bot.axvline(t_pred, color="orange", lw=1.5, ls=":",
                                   label=f"predicted t={t_pred:.1e}")

    ax_bot.set_xlabel("time  t")
    ax_bot.set_ylabel("IPR_trib / IPR_fib")
    ax_bot.set_title("Ratio trib/fib vs time  (cross-over = ratio hits 1.0)")
    ax_bot.grid(True, alpha=0.3)
    ax_bot.legend(fontsize=8, loc="best")
    ax_bot.set_ylim(bottom=0)

    plt.tight_layout()
    plt.savefig(fig_out, dpi=140)
    plt.close()
    print(f"  -> {fig_out}")


def alpha_table_and_plot(
    data: dict[tuple[str, float], dict[str, np.ndarray]],
    csv_out: str = "spreading_exponents.csv",
    fig_out: str = "fig_long_alpha_fit.png",
) -> None:
    """
    Fit alpha for each (chain, lambda) with lambda > 0; write CSV; plot
    fits overlaid on the IPR(t) curves for visual inspection.
    """
    chains = sorted({ch for ch, _ in data.keys()})
    lambdas = sorted({lam for _, lam in data.keys() if lam > 0.0})

    rows: list[tuple[str, float, float, float, int]] = []
    for ch in chains:
        for lam in lambdas:
            key = (ch, lam)
            if key not in data:
                continue
            t = data[key]["time"]
            ipr = data[key]["IPR"]
            alpha, log_A, n = fit_alpha(t, ipr)
            rows.append((ch, lam, alpha, log_A, n))

    with open(csv_out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chain", "lambda", "alpha", "log_A", "n_points"])
        for r in rows:
            w.writerow([r[0], f"{r[1]:.4f}", f"{r[2]:.6f}", f"{r[3]:.6f}", r[4]])
    print(f"  -> {csv_out}")

    print("\n[5] Spreading-exponent fits  IPR(t) ~ t^(-alpha) on late tail")
    print(f"    {'chain':>12} {'lambda':>8} {'alpha':>10} {'pts':>6}")
    print("    " + "-" * 38)
    for ch, lam, alpha, _, n in rows:
        print(f"    {ch:>12} {lam:>8.2f} {alpha:>10.4f} {n:>6d}")

    # Visual: log-log curves with the fitted slope drawn over the late tail.
    fig, axes = plt.subplots(1, len(chains), figsize=(6 * len(chains), 5), dpi=140,
                             squeeze=False)
    cmap = plt.cm.plasma(np.linspace(0.05, 0.95, max(len(lambdas), 1)))
    for ax, ch in zip(axes[0], chains):
        for i, lam in enumerate(lambdas):
            key = (ch, lam)
            if key not in data:
                continue
            t = data[key]["time"]
            ipr = data[key]["IPR"]
            mask = (t > 0) & (ipr > 0)
            ax.loglog(t[mask], ipr[mask], color=cmap[i], lw=1.4,
                      label=f"lambda={lam:.1f}")
            alpha, log_A, n = fit_alpha(t, ipr)
            if np.isfinite(alpha) and n >= 4:
                t_tail = t[mask][-n:]
                fit_line = np.exp(log_A) * t_tail ** (-alpha)
                ax.loglog(t_tail, fit_line, ":", color=cmap[i], lw=2.2)
        ax.set_xlabel("time  t")
        ax.set_ylabel("IPR(t)")
        ax.set_title(f"{ch}  (dotted = late-time t^(-alpha) fit)")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8, loc="best")
    plt.tight_layout()
    plt.savefig(fig_out, dpi=140)
    plt.close()
    print(f"  -> {fig_out}")


# ---------------------------------------------------------------------------
# 6. Entry point
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Analyze long-time DNLS output (ipr_vs_time.csv)."
    )
    ap.add_argument("--csv", default="ipr_vs_time.csv",
                    help="input CSV path (default: ipr_vs_time.csv)")
    ap.add_argument("--no-plots", action="store_true",
                    help="skip plotting (sanity report + alpha CSV only)")
    ap.add_argument(
        "--report", action="store_true",
        help=(
            "Brief 4 verification mode: print verbatim sections [1]-[5] "
            "for the T=10^6, N=1000, λ=1.5 run and generate figure."
        ),
    )
    ap.add_argument(
        "--csv-t5", default=None,
        help="T=10^5 baseline CSV for --report mode (default: none; falls back to in-run snapshot)",
    )
    ap.add_argument(
        "--fig", default="figures/T1e6_lambda1p5_N1000.png",
        help="output figure path for --report mode",
    )
    args = ap.parse_args()

    data = load_csv(args.csv)
    if not data:
        print(f"ERROR: no rows loaded from {args.csv}", file=sys.stderr)
        return 1

    if args.report:
        data5 = load_csv(args.csv_t5) if args.csv_t5 else None
        res = verification_report(data, data5=data5, lam=1.5)
        if not args.no_plots:
            print("\nGenerating verification figure ...")
            plot_verification(
                data,
                fits=res["fits"],
                crossover_t=res["crossover_t"],
                ipr_t6=res["ipr_t6"],
                lam=1.5,
                fig_out=args.fig,
            )
        print("\nDone.")
        return 0

    sanity_report(data)

    if not args.no_plots:
        print("\nGenerating figures ...")
        plot_ipr_vs_t(data)
        plot_lambda0_check(data)

    alpha_table_and_plot(data)
    print("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
