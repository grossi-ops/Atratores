#!/usr/bin/env python3
"""
dnls_lambda1p5_sweep.py
=======================
Targeted λ=1.5 sweep at N∈{500,1000,2000}, T=10⁴.

Directly tests the ΔIPR=−57%/−5% claim from the paper at the published
nonlinearity value.  The existing grid (λ∈{0,1,2,4,8,10}) skips λ=1.5.

Outputs a single long-format CSV ``ipr_lambda1p5_T1e4.csv`` with columns::

    time, lambda, N, chain, IPR, norm

The extra ``N`` column allows a single file to carry all three system sizes.
This file is consumed by ``fss_analyze.py --lambda-sweep`` to incorporate
λ=1.5 into the FSS figures and tables.

Estimated runtime: ~6 minutes on a single core (2 chains × 3 N values).

Usage
-----
Run the full sweep::

    python3 dnls_lambda1p5_sweep.py

Custom N values or additional lambdas::

    python3 dnls_lambda1p5_sweep.py --n-values 500 1000 2000 4000 --lambdas 1.5 2.5

Then fold results into the FSS analysis::

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
import time as _time

import numpy as np

from dnls_long_time import (
    build_hamiltonian,
    evolve_dnls,
    fibonacci_word,
    mid_gap_state,
    tribonacci_word,
)


# ---------------------------------------------------------------------------
# Sweep constants
# ---------------------------------------------------------------------------

N_VALUES: list[int] = [500, 1000, 2000]
LAMBDAS: list[float] = [1.5]
T_END: float = 10_000.0
N_CHECKPOINTS: int = 300
NORM_TOL: float = 1e-5
RTOL: float = 1e-8
ATOL: float = 1e-10
OUT_CSV: str = "ipr_lambda1p5_T1e4.csv"

CHAINS: dict[str, object] = {
    "fibonacci": fibonacci_word,
    "tribonacci": tribonacci_word,
}


# ---------------------------------------------------------------------------
# Sweep
# ---------------------------------------------------------------------------

def run_sweep(
    n_values: list[int] = N_VALUES,
    lambdas: list[float] = LAMBDAS,
    t_end: float = T_END,
    n_checkpoints: int = N_CHECKPOINTS,
    norm_tol: float = NORM_TOL,
    out_csv: str = OUT_CSV,
    verbose: bool = True,
) -> None:
    """
    Sweep over n_values × CHAINS × lambdas, integrate DNLS to t_end,
    and write a combined long-format CSV.

    Output CSV columns
    ------------------
    time    - checkpoint time (log-spaced)
    lambda  - nonlinearity strength
    N       - chain length
    chain   - "fibonacci" or "tribonacci"
    IPR     - inverse participation ratio at that time
    norm    - L2-norm at that time (should remain ~ 1.0)
    """
    rows: list[dict] = []
    n_runs = len(n_values) * len(CHAINS) * len(lambdas)
    run_idx = 0

    for n in n_values:
        for chain_name, word_fn in CHAINS.items():
            word = word_fn(n)  # type: ignore[operator]
            H, hoppings = build_hamiltonian(word, n)
            psi0, E0 = mid_gap_state(H)

            if verbose:
                print(f"\nN={n}  chain={chain_name}  E0={E0:.6f}")

            for lam in lambdas:
                run_idx += 1
                t_wall = _time.perf_counter()

                if verbose:
                    print(
                        f"  [{run_idx}/{n_runs}] lambda={lam}  T={t_end:.0f} ...",
                        end=" ",
                        flush=True,
                    )

                t_arr, ipr_arr, norm_arr, norm_ok = evolve_dnls(
                    psi0,
                    lam,
                    hoppings,
                    t_end=t_end,
                    n_checkpoints=n_checkpoints,
                    norm_tol=norm_tol,
                    rtol=RTOL,
                    atol=ATOL,
                )

                elapsed = _time.perf_counter() - t_wall
                flag = "" if norm_ok else "  *** NORM LEAK ***"
                if verbose:
                    print(
                        f"done in {elapsed:.1f}s  "
                        f"IPR_final={ipr_arr[-1]:.6f}{flag}"
                    )

                for t_k, ipr_k, norm_k in zip(t_arr, ipr_arr, norm_arr):
                    rows.append(
                        {
                            "time": t_k,
                            "lambda": lam,
                            "N": n,
                            "chain": chain_name,
                            "IPR": ipr_k,
                            "norm": norm_k,
                        }
                    )

    with open(out_csv, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["time", "lambda", "N", "chain", "IPR", "norm"]
        )
        writer.writeheader()
        writer.writerows(rows)

    if verbose:
        print(f"\nWrote {len(rows)} rows -> {out_csv}")
        _print_summary(rows, lambdas, n_values)


# ---------------------------------------------------------------------------
# Summary helpers
# ---------------------------------------------------------------------------

def _fit_late_alpha(
    t: np.ndarray,
    ipr: np.ndarray,
    late_frac: float = 0.3,
) -> float:
    """
    Fit IPR(t) ~ t^(-alpha) on the last `late_frac` of the log(t) range.
    Returns alpha (NaN on failure).
    """
    mask = (t > 0) & (ipr > 0)
    t_m = t[mask]
    ipr_m = ipr[mask]
    if len(t_m) < 6:
        return float("nan")
    log_t = np.log(t_m)
    cutoff = log_t[0] + (1.0 - late_frac) * (log_t[-1] - log_t[0])
    sel = log_t >= cutoff
    if sel.sum() < 4:
        return float("nan")
    slope, _ = np.polyfit(log_t[sel], np.log(ipr_m[sel]), 1)
    return float(-slope)


def _print_summary(
    rows: list[dict],
    lambdas: list[float],
    n_values: list[int],
) -> None:
    """
    Print ΔIPR and Δα tables at T_end for each (lambda, N).
    """
    from collections import defaultdict

    # group time-series by (lam, N, chain)
    series: dict[tuple[float, int, str], list[tuple[float, float]]] = defaultdict(list)
    for r in rows:
        key = (float(r["lambda"]), int(r["N"]), r["chain"])
        series[key].append((float(r["time"]), float(r["IPR"])))

    # sort each series by time
    for key in series:
        series[key].sort()

    print()
    print("=" * 70)
    print("SUMMARY  at T=T_end, λ=1.5")
    print()
    print(f"  {'N':>6}  {'IPR_fib':>12}  {'IPR_trib':>12}  {'R=trib/fib':>12}  {'ΔIPR%':>10}")
    print("  " + "-" * 58)
    for lam in sorted(lambdas):
        for n in sorted(n_values):
            f_pts = series.get((lam, n, "fibonacci"), [])
            t_pts = series.get((lam, n, "tribonacci"), [])
            if not f_pts or not t_pts:
                continue
            f_ipr = f_pts[-1][1]
            t_ipr = t_pts[-1][1]
            r = t_ipr / f_ipr
            delta_pct = (t_ipr - f_ipr) / f_ipr * 100.0
            direction = "trib>fib" if r > 1.0 else "trib<fib"
            print(
                f"  {n:>6d}  {f_ipr:>12.7f}  {t_ipr:>12.7f}  "
                f"{r:>12.4f}  {delta_pct:>+10.1f}%  ({direction})"
            )

    print()
    print(f"  {'N':>6}  {'α_fib':>9}  {'α_trib':>9}  {'Δα':>10}  {'Δα/α_fib':>10}")
    print("  " + "-" * 48)
    for lam in sorted(lambdas):
        for n in sorted(n_values):
            f_pts = series.get((lam, n, "fibonacci"), [])
            t_pts = series.get((lam, n, "tribonacci"), [])
            if not f_pts or not t_pts:
                continue
            t_f = np.array([p[0] for p in f_pts])
            ipr_f = np.array([p[1] for p in f_pts])
            t_t = np.array([p[0] for p in t_pts])
            ipr_t = np.array([p[1] for p in t_pts])
            alpha_f = _fit_late_alpha(t_f, ipr_f)
            alpha_t = _fit_late_alpha(t_t, ipr_t)
            if np.isnan(alpha_f) or np.isnan(alpha_t):
                print(f"  {n:>6d}  (fit failed)")
                continue
            da = alpha_t - alpha_f
            da_rel = da / alpha_f * 100.0 if alpha_f != 0 else float("nan")
            direction = "fib spreads faster" if da < 0 else "trib spreads faster"
            print(
                f"  {n:>6d}  {alpha_f:>9.4f}  {alpha_t:>9.4f}  "
                f"{da:>+10.4f}  {da_rel:>+10.1f}%  ({direction})"
            )
    print("=" * 70)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "λ=1.5 sweep at N∈{500,1000,2000}, T=10⁴. "
            "Outputs ipr_lambda1p5_T1e4.csv for use with fss_analyze.py."
        )
    )
    ap.add_argument(
        "--n-values",
        type=int,
        nargs="+",
        default=N_VALUES,
        help=f"chain lengths (default: {N_VALUES})",
    )
    ap.add_argument(
        "--lambdas",
        type=float,
        nargs="+",
        default=LAMBDAS,
        help=f"nonlinearity values (default: {LAMBDAS})",
    )
    ap.add_argument(
        "-T", "--t-end",
        type=float,
        default=T_END,
        help=f"final integration time (default: {T_END})",
    )
    ap.add_argument(
        "--checkpoints",
        type=int,
        default=N_CHECKPOINTS,
        help=f"log-spaced checkpoints in (1, T] (default: {N_CHECKPOINTS})",
    )
    ap.add_argument(
        "--out",
        default=OUT_CSV,
        help=f"output CSV path (default: {OUT_CSV})",
    )
    ap.add_argument(
        "--quiet",
        action="store_true",
        help="suppress progress output",
    )
    args = ap.parse_args()

    run_sweep(
        n_values=args.n_values,
        lambdas=args.lambdas,
        t_end=args.t_end,
        n_checkpoints=args.checkpoints,
        out_csv=args.out,
        verbose=not args.quiet,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
