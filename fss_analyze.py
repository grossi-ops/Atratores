#!/usr/bin/env python3
"""
fss_analyze.py
==============
Finite-size scaling (FSS) analysis companion for the DNLS paper.

Reads a merged long-format CSV produced by running `dnls_long_time.py` at
multiple system sizes N and merging the outputs.  The expected CSV format has
columns:

    N, time, lambda, chain, IPR, norm

Modes
-----
--lambda-sweep FILE
    Analyse the FSS of a single nonlinearity λ across system sizes N.
    Typical use: the T=10⁴ data at λ=1.5 (ipr_lambda1p5_T1e4.csv).

    Produces:
      1. Console report — D₂, spreading exponents α(N), IPR ratios
      2. fig_fss_ipr_vs_t.png       — IPR(t) for both chains at all N
      3. fig_fss_d2_fit.png         — log IPR(t=0) vs log N with D₂ slope
      4. fig_fss_alpha_vs_N.png     — spreading exponent α vs N
      5. fig_fss_ratio_vs_N.png     — trib/fib IPR ratio vs N at T_final

Companion to
------------
  "Differential Nonlinear Robustness of Critical States in Fibonacci and
   Tribonacci Substitution Chains"
  Pablo Nogueira Grossi, G6 LLC (2026)
  DOI: 10.5281/zenodo.20026943

Usage
-----
    python3 fss_analyze.py --lambda-sweep ipr_lambda1p5_T1e4.csv

Author
------
    Pablo Nogueira Grossi  |  ORCID: 0009-0000-6496-2186
    G6 LLC, Newark NJ  |  pablogrossi@hotmail.com

License: MIT
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# 1. CSV loader
# ---------------------------------------------------------------------------

def load_fss_csv(path: str) -> dict[tuple[int, str], dict[str, np.ndarray]]:
    """
    Load the merged FSS CSV (columns: N, time, lambda, chain, IPR, norm).

    Returns {(N, chain): {"time": ..., "IPR": ..., "norm": ..., "lambda": float}}.
    """
    rows_by_key: dict[tuple[int, str], list] = defaultdict(list)
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            key = (int(r["N"]), r["chain"])
            rows_by_key[key].append(
                (float(r["time"]), float(r["IPR"]), float(r["norm"]),
                 float(r["lambda"]))
            )

    out: dict[tuple[int, str], dict] = {}
    for key, rows in rows_by_key.items():
        rows.sort(key=lambda x: x[0])
        out[key] = {
            "time":   np.array([r[0] for r in rows]),
            "IPR":    np.array([r[1] for r in rows]),
            "norm":   np.array([r[2] for r in rows]),
            "lambda": rows[0][3],
        }
    return out


# ---------------------------------------------------------------------------
# 2. Spreading-exponent fit
# ---------------------------------------------------------------------------

def fit_alpha(t: np.ndarray, ipr: np.ndarray, late_frac: float = 0.3
              ) -> tuple[float, float, int]:
    """
    Fit IPR(t) ~ A * t^{-alpha} on the late `late_frac` fraction of log-time.

    Returns (alpha, log_A, n_points).  alpha > 0 means spreading; alpha ~ 0
    indicates self-trapping.
    """
    mask = (t > 0) & (ipr > 0)
    t, ipr = t[mask], ipr[mask]
    if len(t) < 4:
        return float("nan"), float("nan"), 0
    log_t = np.log(t)
    log_ipr = np.log(ipr)
    cutoff = log_t[0] + (1.0 - late_frac) * (log_t[-1] - log_t[0])
    sel = log_t >= cutoff
    if sel.sum() < 3:
        return float("nan"), float("nan"), int(sel.sum())
    slope, intercept = np.polyfit(log_t[sel], log_ipr[sel], 1)
    return float(-slope), float(intercept), int(sel.sum())


# ---------------------------------------------------------------------------
# 3. Console report
# ---------------------------------------------------------------------------

def lambda_sweep_report(data: dict[tuple[int, str], dict]) -> dict:
    """
    Print a structured FSS report and return the computed quantities.
    """
    Ns     = sorted({N for N, _ in data})
    chains = sorted({ch for _, ch in data})

    lam_set = {v["lambda"] for v in data.values()}
    lam = next(iter(lam_set))
    if len(lam_set) > 1:
        print(f"WARNING: multiple lambda values found {lam_set}; using {lam}")

    print("=" * 72)
    print(f"FSS ANALYSIS  —  lambda = {lam},  N ∈ {Ns}")
    print("=" * 72)

    # --- norm conservation ---
    print("\n[0] Norm conservation")
    print(f"    {'N':>6} {'chain':>12} {'max|norm-1|':>14}")
    print("    " + "-" * 34)
    for N in Ns:
        for ch in chains:
            key = (N, ch)
            if key not in data:
                continue
            drift = float(np.max(np.abs(data[key]["norm"] - 1.0)))
            tag = "" if drift < 1e-4 else "  <-- check"
            print(f"    {N:>6} {ch:>12} {drift:>14.2e}{tag}")

    # --- linear-limit D2 (from t=0 IPR scaling) ---
    print("\n[1] Multifractal dimension D₂  (IPR(t=0) ~ N^{-D₂})")
    print("    Uses the t=0 eigenstate IPR — uncontaminated by nonlinear dynamics.")
    print()
    ipr0_by_chain: dict[str, dict[int, float]] = {ch: {} for ch in chains}
    for N in Ns:
        for ch in chains:
            key = (N, ch)
            if key in data:
                ipr0_by_chain[ch][N] = float(data[key]["IPR"][0])

    d2_by_chain: dict[str, float] = {}
    print(f"    {'chain':>12} {'D2':>8}  " + "  ".join(f"IPR(N={n})" for n in Ns))
    print("    " + "-" * (24 + 14 * len(Ns)))
    for ch in chains:
        d = ipr0_by_chain[ch]
        Ns_avail = sorted(d)
        if len(Ns_avail) >= 2:
            log_N = np.log(np.array(Ns_avail, float))
            log_ipr = np.log(np.array([d[n] for n in Ns_avail]))
            slope, _ = np.polyfit(log_N, log_ipr, 1)
            d2 = float(-slope)
        else:
            d2 = float("nan")
        d2_by_chain[ch] = d2
        ipr_str = "  ".join(f"{d.get(n, float('nan')):>12.6f}" for n in Ns)
        print(f"    {ch:>12} {d2:>8.4f}  {ipr_str}")

    # --- spreading exponents ---
    print("\n[2] Spreading exponents α  (IPR(t) ~ t^{-α}, late 30% of log-time)")
    print()
    print(f"    {'N':>6} {'chain':>12} {'alpha':>10} {'n_pts':>6}")
    print("    " + "-" * 38)
    alpha_table: dict[tuple[int, str], float] = {}
    for N in Ns:
        for ch in chains:
            key = (N, ch)
            if key not in data:
                continue
            alpha, _, n = fit_alpha(data[key]["time"], data[key]["IPR"])
            alpha_table[(N, ch)] = alpha
            print(f"    {N:>6} {ch:>12} {alpha:>10.5f} {n:>6d}")

    # --- final IPR and trib/fib ratio ---
    print(f"\n[3] IPR at T_final  (trib/fib ratio)")
    print()
    print(f"    {'N':>6}  " + "  ".join(f"{'IPR_'+ch[:4]:>14}" for ch in chains)
          + "  ratio(t/f)")
    print("    " + "-" * (10 + 16 * len(chains) + 12))
    ratio_table: dict[int, float] = {}
    for N in Ns:
        ipr_final: dict[str, float] = {}
        for ch in chains:
            key = (N, ch)
            if key in data:
                ipr_final[ch] = float(data[key]["IPR"][-1])
        line = f"    {N:>6}  " + "  ".join(
            f"{ipr_final.get(ch, float('nan')):>14.6f}" for ch in chains
        )
        if "fibonacci" in ipr_final and "tribonacci" in ipr_final:
            ratio = ipr_final["tribonacci"] / ipr_final["fibonacci"]
            ratio_table[N] = ratio
            line += f"  {ratio:>10.4f}×"
        print(line)

    print()
    print("=" * 72)
    return {
        "Ns": Ns, "chains": chains, "lam": lam,
        "ipr0_by_chain": ipr0_by_chain,
        "d2_by_chain": d2_by_chain,
        "alpha_table": alpha_table,
        "ratio_table": ratio_table,
    }


# ---------------------------------------------------------------------------
# 4. Figures
# ---------------------------------------------------------------------------

CHAIN_STYLE = {"fibonacci": ("--", 1.3), "tribonacci": ("-", 1.8)}
CHAIN_COLOR = {"fibonacci": "#2980b9", "tribonacci": "#e67e22"}


def fig_ipr_vs_t(data: dict, results: dict, out: str = "fig_fss_ipr_vs_t.png") -> None:
    """IPR(t) for both chains at all N on a shared log-log axes."""
    Ns, chains, lam = results["Ns"], results["chains"], results["lam"]
    cmap = plt.cm.viridis(np.linspace(0.1, 0.9, len(Ns)))

    fig, ax = plt.subplots(figsize=(8, 6), dpi=140)
    for N_idx, N in enumerate(Ns):
        for ch in chains:
            key = (N, ch)
            if key not in data:
                continue
            t = data[key]["time"]
            ipr = data[key]["IPR"]
            mask = t > 0
            ls, lw = CHAIN_STYLE.get(ch, ("-", 1.5))
            ax.loglog(t[mask], ipr[mask], ls, color=cmap[N_idx], lw=lw,
                      label=f"{ch[:4]}  N={N}")

    ax.set_xlabel("time  t")
    ax.set_ylabel("IPR(t)")
    ax.set_title(f"FSS: IPR(t) at λ={lam}  (dashed=Fibonacci, solid=Tribonacci)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(ncol=2, fontsize=9, loc="best")
    plt.tight_layout()
    plt.savefig(out, dpi=140)
    plt.close()
    print(f"  -> {out}")


def fig_d2_fit(results: dict, out: str = "fig_fss_d2_fit.png") -> None:
    """log IPR(t=0) vs log N with fitted D₂ slope for each chain."""
    Ns = results["Ns"]
    ipr0 = results["ipr0_by_chain"]
    d2   = results["d2_by_chain"]
    lam  = results["lam"]

    fig, ax = plt.subplots(figsize=(6, 5), dpi=140)
    log_N_range = np.linspace(np.log(min(Ns)) * 0.95, np.log(max(Ns)) * 1.05, 50)

    for ch in results["chains"]:
        d = ipr0[ch]
        Ns_avail = sorted(d)
        if not Ns_avail:
            continue
        col = CHAIN_COLOR.get(ch, "gray")
        ax.scatter([np.log(n) for n in Ns_avail],
                   [np.log(d[n]) for n in Ns_avail],
                   color=col, zorder=5, s=60, label=f"{ch}  D₂={d2[ch]:.4f}")
        if np.isfinite(d2[ch]) and len(Ns_avail) >= 2:
            log_N_a = np.array([np.log(n) for n in Ns_avail])
            log_ipr_a = np.array([np.log(d[n]) for n in Ns_avail])
            _, intercept = np.polyfit(log_N_a, log_ipr_a, 1)
            fit_y = intercept - d2[ch] * log_N_range
            ax.plot(log_N_range, fit_y, color=col, lw=1.5, alpha=0.6)

    ax.set_xlabel("log N")
    ax.set_ylabel("log IPR(t=0)")
    ax.set_title(f"Multifractal dimension D₂  (λ={lam}, linear limit)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    # annotate x-axis with actual N values
    ax.set_xticks([np.log(n) for n in Ns])
    ax.set_xticklabels([str(n) for n in Ns])
    plt.tight_layout()
    plt.savefig(out, dpi=140)
    plt.close()
    print(f"  -> {out}")


def fig_alpha_vs_N(results: dict, out: str = "fig_fss_alpha_vs_N.png") -> None:
    """Spreading exponent α vs N for both chains."""
    Ns, chains, lam = results["Ns"], results["chains"], results["lam"]
    alpha_table = results["alpha_table"]

    fig, ax = plt.subplots(figsize=(6, 5), dpi=140)
    for ch in chains:
        alphas = [alpha_table.get((N, ch), float("nan")) for N in Ns]
        col = CHAIN_COLOR.get(ch, "gray")
        ls, lw = CHAIN_STYLE.get(ch, ("-", 1.5))
        ax.plot(Ns, alphas, ls + "o", color=col, lw=lw, ms=7,
                label=f"{ch}  (α → 0 means self-trapping)")

    ax.axhline(0, color="gray", lw=0.8, alpha=0.5, ls=":")
    ax.set_xlabel("N  (chain length)")
    ax.set_ylabel("spreading exponent  α")
    ax.set_title(f"α(N) at λ={lam}  — thermodynamic-limit signatures")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    plt.tight_layout()
    plt.savefig(out, dpi=140)
    plt.close()
    print(f"  -> {out}")


def fig_ratio_vs_N(results: dict, out: str = "fig_fss_ratio_vs_N.png") -> None:
    """Trib/Fib IPR ratio at T_final vs N."""
    Ns, lam = results["Ns"], results["lam"]
    ratio_table = results["ratio_table"]

    Ns_avail = sorted(ratio_table)
    ratios = [ratio_table[n] for n in Ns_avail]

    fig, ax = plt.subplots(figsize=(6, 5), dpi=140)
    ax.plot(Ns_avail, ratios, "-o", color="#8e44ad", lw=2, ms=8)
    ax.axhline(1.0, color="gray", lw=0.8, ls=":", alpha=0.6)

    for N, r in zip(Ns_avail, ratios):
        ax.annotate(f"{r:.2f}×", (N, r), textcoords="offset points",
                    xytext=(4, 6), fontsize=9)

    ax.set_xlabel("N  (chain length)")
    ax.set_ylabel("IPR_trib / IPR_fib  at T_final")
    ax.set_title(f"Differential robustness ratio vs N  (λ={lam})")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out, dpi=140)
    plt.close()
    print(f"  -> {out}")


# ---------------------------------------------------------------------------
# 5. Entry point
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description="FSS analysis for DNLS substitution-chain paper."
    )
    ap.add_argument(
        "--lambda-sweep", metavar="FILE",
        help="merged FSS CSV (columns: N,time,lambda,chain,IPR,norm)"
    )
    ap.add_argument("--no-plots", action="store_true",
                    help="skip figure generation")
    ap.add_argument("--late-frac", type=float, default=0.30,
                    help="fraction of late log-time window for alpha fit (default: 0.30)")
    args = ap.parse_args()

    if args.lambda_sweep is None:
        ap.print_help()
        return 1

    data = load_fss_csv(args.lambda_sweep)
    if not data:
        print(f"ERROR: no rows loaded from {args.lambda_sweep}", file=sys.stderr)
        return 1

    results = lambda_sweep_report(data)

    if not args.no_plots:
        print("\nGenerating figures ...")
        fig_ipr_vs_t(data, results)
        fig_d2_fit(results)
        fig_alpha_vs_N(results)
        fig_ratio_vs_N(results)

    print("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
