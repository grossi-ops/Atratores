#!/usr/bin/env python3
"""
fss_analyze.py
==============
Finite-size scaling (FSS) analysis companion for the DNLS Fibonacci/Tribonacci paper.

Companion to:
  "Differential Nonlinear Robustness of Critical States in Fibonacci and
   Tribonacci Substitution Chains"
  Pablo Nogueira Grossi, G6 LLC (2026)
  DOI: 10.5281/zenodo.20026943

Modes
-----
  --csv <path> --report
      Load a long-format IPR CSV (columns: N, time, lambda, chain, IPR, norm)
      and print a structured verification report:
        - N values and chain labels present
        - Final-time IPR values (no NaN, no norm violations)
        - Norm-conservation summary
        - Differential ratio IPR_fib / IPR_trib at final time

  --lambda-sweep [--csv <path>]
      Read the full FSS dataset (default: data/ipr_fss_T1e4.csv) and print
      the complete FSS table: IPR_fib/IPR_trib at final time, all (lambda, N).

Usage
-----
    python3 fss_analyze.py --csv data/ipr_lambda1p5_T1e4.csv --report
    python3 fss_analyze.py --lambda-sweep
    python3 fss_analyze.py --lambda-sweep --csv data/ipr_fss_T1e4.csv

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


# ---------------------------------------------------------------------------
# 1. CSV loader
# ---------------------------------------------------------------------------

def load_csv(path: str) -> dict[tuple[int, float, str], dict[str, np.ndarray]]:
    """
    Load a long-format IPR CSV into a nested dict.

    Required columns: N (or n_sites), time, lambda, chain, IPR, norm.
    Returns {(N, lambda, chain): {"time": arr, "IPR": arr, "norm": arr}}.
    """
    rows_by_key: dict[tuple[int, float, str], list[tuple[float, float, float]]] = (
        defaultdict(list)
    )
    with open(path, "r", newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames or []
        n_col = "N" if "N" in fieldnames else "n_sites"
        for r in reader:
            key = (int(float(r[n_col])), float(r["lambda"]), r["chain"])
            rows_by_key[key].append(
                (float(r["time"]), float(r["IPR"]), float(r["norm"]))
            )

    out: dict[tuple[int, float, str], dict[str, np.ndarray]] = {}
    for key, rows in rows_by_key.items():
        rows.sort(key=lambda x: x[0])
        out[key] = {
            "time": np.array([r[0] for r in rows]),
            "IPR":  np.array([r[1] for r in rows]),
            "norm": np.array([r[2] for r in rows]),
        }
    return out


# ---------------------------------------------------------------------------
# 2. Verification report  (--csv --report)
# ---------------------------------------------------------------------------

NORM_WARN = 5e-5   # flag if max |norm - 1| exceeds this


def verification_report(
    data: dict[tuple[int, float, str], dict[str, np.ndarray]],
    norm_warn: float = NORM_WARN,
) -> dict:
    """
    Print a structured verification report and return a summary dict.

    Checks:
      1. Which (N, lambda, chain) tuples are present
      2. NaN / Inf in final-time IPR
      3. Norm conservation: max |norm - 1| per run
      4. Differential ratio IPR_fib / IPR_trib at final time
    """
    Ns      = sorted({k[0] for k in data})
    lambdas = sorted({k[1] for k in data})
    chains  = sorted({k[2] for k in data})

    print("=" * 72)
    print("FSS VERIFICATION REPORT")
    print("=" * 72)

    # ── Section 1: coverage ────────────────────────────────────────────────
    print(f"\n[1] Coverage")
    print(f"    N values  : {Ns}")
    print(f"    λ values  : {lambdas}")
    print(f"    chains    : {chains}")
    missing = []
    for N in Ns:
        for lam in lambdas:
            for ch in chains:
                if (N, lam, ch) not in data:
                    missing.append((N, lam, ch))
    if missing:
        print(f"    MISSING   : {missing}")
    else:
        print(f"    No missing (N, λ, chain) combinations.")

    # ── Section 2: final-time IPR sanity ───────────────────────────────────
    print(f"\n[2] Final-time IPR sanity (NaN / Inf check)")
    bad_ipr = []
    for key, d in data.items():
        v = d["IPR"][-1]
        if not np.isfinite(v):
            bad_ipr.append((key, v))
    if bad_ipr:
        for key, v in bad_ipr:
            print(f"    BAD IPR: N={key[0]}  λ={key[1]}  chain={key[2]}  IPR={v}")
    else:
        print(f"    All final-time IPR values finite: OK")

    # ── Section 3: norm conservation ───────────────────────────────────────
    print(f"\n[3] Norm conservation  (flag threshold: {norm_warn:.1e})")
    print(f"    {'N':>5}  {'λ':>5}  {'chain':>12}  {'max|norm-1|':>14}  {'pts':>5}")
    print("    " + "-" * 50)
    norm_summary: list[tuple[int, float, str, float]] = []
    worst = 0.0
    for (N, lam, ch), d in sorted(data.items()):
        drift = float(np.max(np.abs(d["norm"] - 1.0)))
        worst = max(worst, drift)
        flag = "  ← FLAG" if drift > norm_warn else ""
        print(f"    {N:>5}  {lam:>5.1f}  {ch:>12}  {drift:>14.2e}  {len(d['time']):>5}{flag}")
        norm_summary.append((N, lam, ch, drift))
    print(f"    Worst drift: {worst:.2e}{'  ← FLAG' if worst > norm_warn else '  OK'}")

    # ── Section 4: differential ratio at final time ────────────────────────
    print(f"\n[4] Differential ratio IPR_fib / IPR_trib at final time T")
    Ns_with_both = [
        N for N in Ns
        if any((N, lam, "fibonacci") in data for lam in lambdas)
        and any((N, lam, "tribonacci") in data for lam in lambdas)
    ]
    for lam in lambdas:
        row = f"    λ={lam:.1f}:"
        for N in Ns:
            fk = (N, lam, "fibonacci")
            tk = (N, lam, "tribonacci")
            if fk in data and tk in data:
                f_ipr = float(data[fk]["IPR"][-1])
                t_ipr = float(data[tk]["IPR"][-1])
                ratio = f_ipr / t_ipr if t_ipr > 0 else float("nan")
                row += f"  N={N}: {ratio:.4f}"
        print(row)

    print("\n" + "=" * 72)
    return {
        "Ns": Ns,
        "lambdas": lambdas,
        "chains": chains,
        "missing": missing,
        "bad_ipr": bad_ipr,
        "norm_worst": worst,
        "norm_summary": norm_summary,
        "norm_ok": worst <= norm_warn,
    }


# ---------------------------------------------------------------------------
# 3. Lambda-sweep FSS table  (--lambda-sweep)
# ---------------------------------------------------------------------------

def lambda_sweep_report(
    data: dict[tuple[int, float, str], dict[str, np.ndarray]],
) -> None:
    """
    Print the complete FSS table: IPR_fib/IPR_trib at final time.

    Format (3 significant figures, fixed-width):

        λ        N=500    N=1000    N=2000
        ─────────────────────────────────
        0.5      0.171     0.136     0.113
        1.0      0.615     0.343     0.159
        1.5      0.658     0.239     0.315   ← non-monotone
        ...
    """
    lambdas = sorted({k[1] for k in data})
    Ns      = sorted({k[0] for k in data})

    header = f"    {'λ':>5}" + "".join(f"  {('N='+str(N)):>8}" for N in Ns)
    print()
    print("FSS TABLE — differential ratio IPR_fib/IPR_trib at final time T")
    print()
    print(header)
    print("    " + "─" * (len(header) - 4))

    for lam in lambdas:
        row = f"    {lam:>5.1f}"
        ratios = []
        for N in Ns:
            fk = (N, lam, "fibonacci")
            tk = (N, lam, "tribonacci")
            if fk in data and tk in data:
                f_ipr = float(data[fk]["IPR"][-1])
                t_ipr = float(data[tk]["IPR"][-1])
                ratio = f_ipr / t_ipr if t_ipr > 0 else float("nan")
                ratios.append(ratio)
                row += f"  {ratio:>8.3f}"
            else:
                ratios.append(float("nan"))
                row += f"  {'N/A':>8}"
        # Flag non-monotone rows
        finite = [r for r in ratios if np.isfinite(r)]
        if len(finite) >= 3:
            diffs = [finite[i+1] - finite[i] for i in range(len(finite) - 1)]
            if not all(d > 0 for d in diffs) and not all(d < 0 for d in diffs):
                row += "  ← non-monotone"
        print(row)


# ---------------------------------------------------------------------------
# 4. Entry point
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description="FSS analysis for DNLS Fibonacci/Tribonacci chains."
    )
    ap.add_argument(
        "--csv",
        default="data/ipr_fss_T1e4.csv",
        help="Input CSV path (default: data/ipr_fss_T1e4.csv)",
    )
    ap.add_argument(
        "--report",
        action="store_true",
        help="Run full verification report on the given --csv file.",
    )
    ap.add_argument(
        "--lambda-sweep",
        action="store_true",
        dest="lambda_sweep",
        help="Print the FSS table (IPR_fib/IPR_trib) for all (λ, N) in --csv.",
    )
    args = ap.parse_args()

    if not args.report and not args.lambda_sweep:
        ap.print_help()
        return 1

    try:
        data = load_csv(args.csv)
    except FileNotFoundError:
        print(f"ERROR: file not found: {args.csv}", file=sys.stderr)
        return 1

    if not data:
        print(f"ERROR: no rows loaded from {args.csv}", file=sys.stderr)
        return 1

    if args.report:
        verification_report(data)

    if args.lambda_sweep:
        lambda_sweep_report(data)

    return 0


if __name__ == "__main__":
    sys.exit(main())
