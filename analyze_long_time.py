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

Companion to:
  "Differential Nonlinear Robustness of Critical States in Fibonacci and
   Tribonacci Substitution Chains"
  Pablo Nogueira Grossi, G6 LLC (2026)
  DOI: 10.5281/zenodo.20026943

Usage
-----
    python3 analyze_long_time.py                          # uses ipr_vs_time.csv
    python3 analyze_long_time.py --csv long_run.csv       # custom input

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
    t_min_fit: float = 1e4,
) -> tuple[float, float, float, float, int, float]:
    """
    Fit IPR(t) ~ A * t^(-alpha) for t >= t_min_fit using OLS with stderr and R².

    Returns
    -------
    alpha       : spreading exponent (positive means spreading)
    alpha_stderr: standard error of alpha from OLS
    log_A       : intercept in log space
    R2          : coefficient of determination in log space
    n_pts       : number of points used
    t_min_used  : actual minimum time used (first checkpoint >= t_min_fit)
    """
    mask = (t >= t_min_fit) & (ipr > 0)
    n = int(mask.sum())
    if n < 4:
        return float("nan"), float("nan"), float("nan"), float("nan"), n, float("nan")

    log_t = np.log(t[mask])
    log_ipr = np.log(ipr[mask])
    t_min_used = float(t[mask][0])

    # OLS: log_ipr = intercept + slope * log_t
    A_mat = np.column_stack([np.ones(n), log_t])
    coeffs, residuals, rank, _ = np.linalg.lstsq(A_mat, log_ipr, rcond=None)
    intercept, slope = float(coeffs[0]), float(coeffs[1])

    # residuals sum-of-squares
    log_ipr_pred = intercept + slope * log_t
    ss_res = float(np.sum((log_ipr - log_ipr_pred) ** 2))
    ss_tot = float(np.sum((log_ipr - log_ipr.mean()) ** 2))
    R2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    # standard error of slope
    if n > 2:
        s2 = ss_res / (n - 2)
        ss_xx = float(np.sum((log_t - log_t.mean()) ** 2))
        slope_se = float(np.sqrt(s2 / ss_xx)) if ss_xx > 0 else float("nan")
    else:
        slope_se = float("nan")

    return float(-slope), slope_se, intercept, R2, n, t_min_used


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


def alpha_fit_table(
    data: dict[tuple[str, float], dict[str, np.ndarray]],
    t_min_fit: float = 1e4,
) -> list[dict]:
    """
    Fit alpha for all (lambda, chain) pairs using t >= t_min_fit window.

    Returns a list of dicts with keys: lambda, chain, alpha, alpha_stderr,
    t_min_fit, R2.
    """
    lambdas = sorted({lam for _, lam in data.keys() if lam > 0.0})
    chains = sorted({ch for ch, _ in data.keys()})
    results = []
    for lam in lambdas:
        for ch in chains:
            key = (ch, lam)
            if key not in data:
                continue
            t = data[key]["time"]
            ipr_arr = data[key]["IPR"]
            alpha, stderr, log_A, R2, n, t_min_used = fit_alpha_full(
                t, ipr_arr, t_min_fit=t_min_fit
            )
            results.append(
                {
                    "lambda": lam,
                    "chain": ch,
                    "alpha": alpha,
                    "alpha_stderr": stderr,
                    "t_min_fit": t_min_used,
                    "R2": R2,
                    "n_pts": n,
                    "log_A": log_A,
                }
            )
    return results


def plot_alpha_N2000(
    data: dict[tuple[str, float], dict[str, np.ndarray]],
    fit_results: list[dict],
    out_path: str = "figures/alpha_N2000_T1e5.png",
    t_min_fit: float = 1e4,
) -> None:
    """
    IPR(t) log-log for all 8 (lambda, chain) pairs with alpha-fit lines
    overlaid in the t > t_min_fit window.  Saved to `out_path`.
    """
    import os
    os.makedirs(os.path.dirname(out_path) if os.path.dirname(out_path) else ".", exist_ok=True)

    lambdas = sorted({lam for _, lam in data.keys() if lam > 0.0})
    chains = sorted({ch for ch, _ in data.keys() if ch in ("fibonacci", "tribonacci")})

    # Build fit lookup
    fit_lut: dict[tuple[float, str], dict] = {
        (r["lambda"], r["chain"]): r for r in fit_results
    }

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), dpi=140, sharey=True)
    cmap = plt.cm.tab10(np.linspace(0, 0.8, max(len(lambdas), 1)))
    ls_map = {"fibonacci": "--", "tribonacci": "-"}
    lw_map = {"fibonacci": 1.4, "tribonacci": 1.9}

    for ax, ch in zip(axes, chains):
        for i, lam in enumerate(lambdas):
            key = (ch, lam)
            if key not in data:
                continue
            t = data[key]["time"]
            ipr_arr = data[key]["IPR"]
            mask = t > 0
            ax.loglog(
                t[mask], ipr_arr[mask],
                ls_map.get(ch, "-"),
                color=cmap[i], lw=lw_map.get(ch, 1.5),
                label=f"λ={lam:.1f}",
            )
            # overlay fit line in the t > t_min_fit window
            fit = fit_lut.get((lam, ch))
            if fit and np.isfinite(fit["alpha"]) and np.isfinite(fit["log_A"]):
                t_fit = t[(t >= t_min_fit) & mask]
                if len(t_fit) >= 2:
                    fit_line = np.exp(fit["log_A"]) * t_fit ** (-fit["alpha"])
                    ax.loglog(
                        t_fit, fit_line,
                        ":", color=cmap[i], lw=2.4,
                    )
                    # annotate alpha value near end of fit line
                    ax.annotate(
                        f"α={fit['alpha']:.3f}",
                        xy=(t_fit[-1], fit_line[-1]),
                        xytext=(4, 2),
                        textcoords="offset points",
                        fontsize=6.5,
                        color=cmap[i],
                    )
        ax.axvline(t_min_fit, color="gray", lw=0.8, ls="--", alpha=0.6,
                   label=f"t_min={t_min_fit:.0e}")
        ax.set_xlabel("time  t")
        ax.set_ylabel("IPR(t)")
        ax.set_title(f"{ch.capitalize()}  N=2000, T=10⁵")
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(fontsize=7, loc="lower left")

    plt.suptitle(
        "IPR(t) log-log — N=2000, T=10⁵ (dotted lines = α-fit for t > 10⁴)",
        fontsize=10,
    )
    plt.tight_layout()
    plt.savefig(out_path, dpi=140)
    plt.close()
    print(f"  -> {out_path}")


# ---------------------------------------------------------------------------
# 5. Entry point
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Analyze long-time DNLS output (ipr_vs_time.csv)."
    )
    ap.add_argument("--csv", default="ipr_vs_time.csv",
                    help="input CSV path (default: ipr_vs_time.csv)")
    ap.add_argument("--no-plots", action="store_true",
                    help="skip plotting (sanity report + alpha CSV only)")
    ap.add_argument("--report", action="store_true",
                    help=(
                        "print the verbatim PR-report sections [1]–[5]: sanity, "
                        "alpha-fit table (t>1e4, with stderr and R²), comparison "
                        "notes, norm conservation summary, and generate "
                        "figures/alpha_N2000_T1e5.png"
                    ))
    ap.add_argument("--t-min-fit", type=float, default=1e4,
                    help="lower time cutoff for alpha regression (default: 1e4)")
    args = ap.parse_args()

    data = load_csv(args.csv)
    if not data:
        print(f"ERROR: no rows loaded from {args.csv}", file=sys.stderr)
        return 1

    report_data = sanity_report(data)

    if args.report:
        # --- [2] α-fit table ---
        t_min_fit = args.t_min_fit
        fit_results = alpha_fit_table(data, t_min_fit=t_min_fit)
        print()
        print("=" * 72)
        print("[2] α-FIT TABLE  (t > {:.0e}, OLS in log-log)".format(t_min_fit))
        print("=" * 72)
        header = (
            f"{'λ':>5}  {'chain':>12}  {'α':>8}  {'α_stderr':>10}  "
            f"{'t_min_fit':>12}  {'R²':>8}"
        )
        print(header)
        print("-" * len(header))
        for r in fit_results:
            alpha_s = f"{r['alpha']:.4f}" if np.isfinite(r["alpha"]) else "   nan"
            se_s = f"{r['alpha_stderr']:.4f}" if np.isfinite(r["alpha_stderr"]) else "      nan"
            tm_s = f"{r['t_min_fit']:.2e}" if np.isfinite(r["t_min_fit"]) else "        nan"
            r2_s = f"{r['R2']:.4f}" if np.isfinite(r["R2"]) else "     nan"
            print(
                f"{r['lambda']:>5.1f}  {r['chain']:>12}  "
                f"{alpha_s:>8}  {se_s:>10}  {tm_s:>12}  {r2_s:>8}"
            )

        # --- [3] Comparison note ---
        print()
        print("=" * 72)
        print("[3] COMPARISON TO T=10⁴ VALUES")
        print("=" * 72)
        print("  α values from ipr_lambda1p5_N1000_T1e5.csv (N=1000, T=10⁵, λ=1.5):")
        print("    fibonacci  λ=1.5 : see previous PR / section8_draft.md")
        print("    tribonacci λ=1.5 : see previous PR / section8_draft.md")
        print()
        print("  New N=2000, T=10⁵ values (this run):")
        for r in fit_results:
            alpha_s = f"{r['alpha']:.4f}" if np.isfinite(r["alpha"]) else "nan"
            print(f"    {r['chain']:>12}  λ={r['lambda']:.1f}  α={alpha_s}")
        fib_15 = next(
            (r for r in fit_results if r["chain"] == "fibonacci" and r["lambda"] == 1.5),
            None,
        )
        trib_15 = next(
            (r for r in fit_results if r["chain"] == "tribonacci" and r["lambda"] == 1.5),
            None,
        )
        if fib_15 and np.isfinite(fib_15["alpha"]):
            in_range = 0.20 <= fib_15["alpha"] <= 0.22
            print(f"\n  α_fib(λ=1.5) = {fib_15['alpha']:.4f}  "
                  f"→ {'IN' if in_range else 'OUTSIDE'} [0.20, 0.22]")
        if trib_15 and fib_15 and np.isfinite(trib_15["alpha"]) and np.isfinite(fib_15["alpha"]):
            decreasing = trib_15["alpha"] < fib_15["alpha"]
            print(f"  α_trib(λ=1.5) = {trib_15['alpha']:.4f}  "
                  f"→ α_trib {'<' if decreasing else '>='} α_fib  "
                  f"(sign-flip {'confirmed' if decreasing else 'NOT confirmed'})")

        # --- [4] Norm conservation ---
        print()
        print("=" * 72)
        print("[4] NORM CONSERVATION")
        print("=" * 72)
        worst = report_data["worst_drift"]
        print(f"  Maximum norm leak across all checkpoints: {worst:.2e}")
        if worst > 5e-5:
            print("  *** EXCEEDS 5e-5 THRESHOLD ***")
        else:
            print("  OK — within 5e-5 threshold.")

        # --- [5] Plot ---
        if not args.no_plots:
            print()
            print("=" * 72)
            print("[5] FIGURE")
            print("=" * 72)
            out_fig = "figures/alpha_N2000_T1e5.png"
            plot_alpha_N2000(data, fit_results, out_path=out_fig, t_min_fit=t_min_fit)

        print()
        print("=" * 72)
        return 0

    if not args.no_plots:
        print("\nGenerating figures ...")
        plot_ipr_vs_t(data)
        plot_lambda0_check(data)

    alpha_table_and_plot(data)
    print("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
