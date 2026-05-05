#!/usr/bin/env python3
"""
fss_analyze.py
==============
Finite-size scaling (FSS) companion for `dnls_long_time.py`.

Reads per-N CSV files produced by running::

    python3 dnls_long_time.py -N <N> -T 10000 --checkpoints 300 --out ipr_N<N>_T1e4.csv

plus the existing N=500 slice extracted from ``ipr_vs_time.csv``.

Expected input files (any subset is accepted; missing sizes are skipped):
    ipr_N200_T1e4.csv
    ipr_N500_T1e4.csv
    ipr_N1000_T1e4.csv
    ipr_N2000_T1e4.csv

Three things looked for
-----------------------
1. **D2 (multifractal dimension)** – linear-limit (λ=0) IPR(N) ~ N^(-D2).
   A critical eigenstate has 0 < D2 < 1; an extended state has D2=1;
   a localised state has D2=0.

2. **Spreading exponent α(N)** – late-time fit IPR(t) ~ t^(-α) for each
   (chain, λ, N).  If α is N-independent the thermodynamic limit is reached;
   a downward drift with N signals finite-size saturation.

3. **Saturation crossover time t_sat(N)** – estimated as the time where
   IPR(t) first falls within a factor of 2 of its final value.  A scaling
   t_sat ~ N^z would reveal the dynamical exponent z.

Outputs
-------
  fss_D2.csv                  – D2 fits per chain
  fss_alpha.csv               – α per (chain, lambda, N)
  fss_tsat.csv                – t_sat per (chain, lambda, N)
  fig_fss_D2.png              – IPR(N) at λ=0 log-log with D2 fit
  fig_fss_alpha_vs_N.png      – α(N) curves per (chain, lambda)
  fig_fss_tsat_vs_N.png       – t_sat(N) curves

Author
------
    Pablo Nogueira Grossi  |  ORCID: 0009-0000-6496-2186
    G6 LLC, Newark NJ  |  GitHub: https://github.com/TOTOGT/AXLE

License: MIT
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# 1. I/O helpers
# ---------------------------------------------------------------------------

def load_csv(path: str) -> dict[tuple[str, float], dict[str, np.ndarray]]:
    """
    Load a long-format CSV (time, lambda, chain, IPR, norm) into
    {(chain, lambda): {time, IPR, norm}}.
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


def slice_at_t(
    data: dict[tuple[str, float], dict[str, np.ndarray]],
    t_max: float,
) -> dict[tuple[str, float], dict[str, np.ndarray]]:
    """Return a view keeping only t <= t_max."""
    out = {}
    for key, arrs in data.items():
        mask = arrs["time"] <= t_max
        out[key] = {k: v[mask] for k, v in arrs.items()}
    return out


def load_all(
    n_values: list[int],
    t_max: float,
) -> dict[int, dict[tuple[str, float], dict[str, np.ndarray]]]:
    """
    Load all available per-N CSVs.  For N=500 the large T=10^5 file is used
    and sliced at t_max.  Returns {N: data_dict}.
    """
    out: dict[int, dict] = {}
    for n in n_values:
        candidate = f"ipr_N{n}_T1e4.csv"
        fallback = "ipr_vs_time.csv" if n == 500 else None
        if os.path.exists(candidate):
            data = load_csv(candidate)
            data = slice_at_t(data, t_max)
            out[n] = data
            print(f"  loaded {candidate}")
        elif fallback and os.path.exists(fallback):
            data = load_csv(fallback)
            data = slice_at_t(data, t_max)
            out[n] = data
            print(f"  loaded {fallback} (sliced at t<={t_max:.0f}) for N={n}")
        else:
            print(f"  [skip] no CSV found for N={n}")
    return out


def load_lambda_sweep(
    path: str,
    t_max: float,
) -> dict[int, dict[tuple[str, float], dict[str, np.ndarray]]]:
    """
    Load a combined multi-N sweep CSV produced by ``dnls_lambda1p5_sweep.py``.

    Expected columns: time, lambda, N, chain, IPR, norm.
    The extra ``N`` column distinguishes this format from the per-N files.

    Returns the same ``{N: {(chain, lambda): arrays}}`` structure used by
    :func:`load_all`, so the two results can be merged with a simple
    ``dict.update``.
    """
    rows_by_key: dict[
        tuple[int, str, float], list[tuple[float, float, float]]
    ] = defaultdict(list)

    with open(path, "r", newline="") as fh:
        reader = csv.DictReader(fh)
        if "N" not in (reader.fieldnames or []):
            raise ValueError(
                f"{path!r} is missing the 'N' column. "
                "Expected a multi-N sweep CSV from dnls_lambda1p5_sweep.py."
            )
        for r in reader:
            t = float(r["time"])
            if t > t_max:
                continue
            n = int(r["N"])
            chain = r["chain"]
            lam = float(r["lambda"])
            rows_by_key[(n, chain, lam)].append(
                (t, float(r["IPR"]), float(r["norm"]))
            )

    out: dict[int, dict[tuple[str, float], dict[str, np.ndarray]]] = defaultdict(dict)
    for (n, chain, lam), rows in rows_by_key.items():
        rows.sort(key=lambda x: x[0])
        out[n][(chain, lam)] = {
            "time": np.array([r[0] for r in rows]),
            "IPR":  np.array([r[1] for r in rows]),
            "norm": np.array([r[2] for r in rows]),
        }

    return dict(out)


# ---------------------------------------------------------------------------
# 2. Fitting helpers
# ---------------------------------------------------------------------------

def fit_power_law(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """
    Fit y ~ A * x^exponent via log-log linear regression.
    Returns (exponent, log_A).
    """
    mask = (x > 0) & (y > 0)
    log_x = np.log(x[mask])
    log_y = np.log(y[mask])
    slope, intercept = np.polyfit(log_x, log_y, 1)
    return float(slope), float(intercept)


def fit_alpha(t: np.ndarray, ipr: np.ndarray, late_frac: float = 0.3
              ) -> tuple[float, float, int]:
    """
    Fit IPR(t) ~ A * t^(-alpha) on the last `late_frac` of the log(t) range.
    Returns (alpha, log_A, n_points).
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


def estimate_tsat(t: np.ndarray, ipr: np.ndarray, factor: float = 2.0) -> float:
    """
    Estimate the saturation crossover time as the first t where
    IPR(t) < factor * IPR_final.  Returns nan if IPR never approaches final.
    """
    mask = t > 0
    t = t[mask]
    ipr = ipr[mask]
    if len(ipr) < 4:
        return float("nan")
    ipr_final = float(ipr[-1])
    threshold = factor * ipr_final
    # walk forward until IPR drops below threshold
    idx = np.argmax(ipr <= threshold)
    if idx == 0 and ipr[0] > threshold:
        return float("nan")   # never reached
    if idx == 0:
        return float(t[0])
    return float(t[idx])


# ---------------------------------------------------------------------------
# 3. Analysis sections
# ---------------------------------------------------------------------------

def analyze_D2(
    all_data: dict[int, dict[tuple[str, float], dict[str, np.ndarray]]],
    chains: list[str],
) -> dict[str, tuple[float, float]]:
    """
    [1] D2 multifractal dimension from IPR(N) at lambda=0 at final time.
    Returns {chain: (D2, log_A)}.
    """
    print("\n" + "=" * 72)
    print("[1] Multifractal dimension D2  --  linear limit (lambda=0)")
    print("    IPR(N) at lambda=0, t=T_max  ~  N^(-D2)")
    print("    D2 = 1 => extended;  D2 = 0 => localised;  0 < D2 < 1 => critical")
    print()
    print(f"    {'chain':>12} {'N':>6} {'IPR(T,lambda=0)':>18}")
    print("    " + "-" * 38)

    result = {}
    for ch in chains:
        ns = sorted(all_data.keys())
        ns_with_data = []
        ipr_vals = []
        for n in ns:
            key = (ch, 0.0)
            if key not in all_data[n]:
                continue
            ipr = float(all_data[n][key]["IPR"][-1])
            print(f"    {ch:>12} {n:>6} {ipr:>18.10f}")
            ns_with_data.append(n)
            ipr_vals.append(ipr)
        if len(ns_with_data) >= 2:
            exp, log_A = fit_power_law(
                np.array(ns_with_data, dtype=float),
                np.array(ipr_vals),
            )
            D2 = float(-exp)
            result[ch] = (D2, log_A)
            print(f"\n    {ch:>12}  =>  D2 = {D2:.4f}  (IPR ~ N^(-{D2:.4f}))\n")
        else:
            print(f"    {ch:>12}  =>  insufficient N values for D2 fit\n")
    return result


def analyze_alpha(
    all_data: dict[int, dict[tuple[str, float], dict[str, np.ndarray]]],
    chains: list[str],
    lambdas: list[float],
) -> list[tuple[str, float, int, float, int]]:
    """
    [2] Spreading exponent alpha(N) for lambda > 0.
    Returns list of (chain, lambda, N, alpha, n_pts).
    """
    print("=" * 72)
    print("[2] Spreading exponent alpha(N)  --  IPR(t) ~ t^(-alpha) late tail")
    print("    N-independence => thermodynamic limit reached.")
    print("    Downward drift with N => finite-size saturation.")
    print()
    print(f"    {'chain':>12} {'lambda':>8} {'N':>6} {'alpha':>10} {'pts':>6}")
    print("    " + "-" * 44)

    rows = []
    for ch in chains:
        for lam in lambdas:
            if lam == 0.0:
                continue
            for n in sorted(all_data.keys()):
                key = (ch, lam)
                if key not in all_data[n]:
                    continue
                t = all_data[n][key]["time"]
                ipr = all_data[n][key]["IPR"]
                alpha, _, n_pts = fit_alpha(t, ipr)
                print(f"    {ch:>12} {lam:>8.2f} {n:>6} {alpha:>10.4f} {n_pts:>6d}")
                rows.append((ch, lam, n, alpha, n_pts))
    print()
    return rows


def analyze_tsat(
    all_data: dict[int, dict[tuple[str, float], dict[str, np.ndarray]]],
    chains: list[str],
    lambdas: list[float],
) -> list[tuple[str, float, int, float]]:
    """
    [3] Saturation crossover time t_sat(N).
    """
    print("=" * 72)
    print("[3] Saturation crossover t_sat(N)  -- first t where IPR < 2 * IPR(T)")
    print("    A scaling t_sat ~ N^z reveals the dynamical exponent z.")
    print()
    print(f"    {'chain':>12} {'lambda':>8} {'N':>6} {'t_sat':>12}")
    print("    " + "-" * 40)

    rows = []
    for ch in chains:
        for lam in lambdas:
            if lam == 0.0:
                continue
            for n in sorted(all_data.keys()):
                key = (ch, lam)
                if key not in all_data[n]:
                    continue
                t = all_data[n][key]["time"]
                ipr = all_data[n][key]["IPR"]
                tsat = estimate_tsat(t, ipr)
                tag = f"{tsat:>12.2f}" if np.isfinite(tsat) else f"{'> T_max':>12}"
                print(f"    {ch:>12} {lam:>8.2f} {n:>6} {tag}")
                rows.append((ch, lam, n, tsat))
    print()
    return rows


# ---------------------------------------------------------------------------
# 4. Plots
# ---------------------------------------------------------------------------

def plot_D2(
    all_data: dict[int, dict[tuple[str, float], dict[str, np.ndarray]]],
    chains: list[str],
    d2_fits: dict[str, tuple[float, float]],
    out_path: str = "fig_fss_D2.png",
) -> None:
    fig, ax = plt.subplots(figsize=(7, 5), dpi=140)
    colors = {"fibonacci": "steelblue", "tribonacci": "darkorange"}
    markers = {"fibonacci": "o", "tribonacci": "s"}

    for ch in chains:
        ns = sorted(all_data.keys())
        ns_ok, ipr_ok = [], []
        for n in ns:
            key = (ch, 0.0)
            if key not in all_data[n]:
                continue
            ns_ok.append(n)
            ipr_ok.append(float(all_data[n][key]["IPR"][-1]))
        if not ns_ok:
            continue
        ax.loglog(ns_ok, ipr_ok, markers.get(ch, "o"), color=colors.get(ch, "gray"),
                  ms=8, label=f"{ch}  data")
        if ch in d2_fits:
            D2, log_A = d2_fits[ch]
            n_arr = np.array([min(ns_ok), max(ns_ok)])
            ax.loglog(n_arr, np.exp(log_A) * n_arr ** (-D2), "--",
                      color=colors.get(ch, "gray"), lw=1.5,
                      label=f"{ch}  fit D2={D2:.3f}")

    ax.set_xlabel("chain length  N")
    ax.set_ylabel("IPR(N, T=T_max, lambda=0)")
    ax.set_title("D2 multifractal dimension  (linear limit, lambda=0)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="best")
    plt.tight_layout()
    plt.savefig(out_path, dpi=140)
    plt.close()
    print(f"  -> {out_path}")


def plot_alpha_vs_N(
    all_data: dict[int, dict[tuple[str, float], dict[str, np.ndarray]]],
    alpha_rows: list[tuple[str, float, int, float, int]],
    chains: list[str],
    lambdas: list[float],
    out_path: str = "fig_fss_alpha_vs_N.png",
) -> None:
    lams_nonzero = [lam for lam in lambdas if lam > 0.0]
    if not lams_nonzero:
        return
    cmap = plt.cm.viridis(np.linspace(0.05, 0.95, len(lams_nonzero)))
    style = {"fibonacci": "--", "tribonacci": "-"}
    marker = {"fibonacci": "^", "tribonacci": "o"}

    fig, ax = plt.subplots(figsize=(8, 5), dpi=140)
    for ch in chains:
        for i, lam in enumerate(lams_nonzero):
            pts = [(n, a) for (c, l, n, a, _) in alpha_rows
                   if c == ch and l == lam and np.isfinite(a)]
            if not pts:
                continue
            pts.sort()
            ns, alphas = zip(*pts)
            ax.plot(ns, alphas, style.get(ch, "-") + marker.get(ch, "o"),
                    color=cmap[i], lw=1.4, ms=6,
                    label=f"{ch[:4]}  lambda={lam:.1f}")

    ax.set_xlabel("chain length  N")
    ax.set_ylabel("spreading exponent  alpha")
    ax.set_title("alpha(N) -- finite-size dependence of spreading exponent")
    ax.grid(True, alpha=0.3)
    ax.legend(ncol=2, fontsize=8, loc="best")
    plt.tight_layout()
    plt.savefig(out_path, dpi=140)
    plt.close()
    print(f"  -> {out_path}")


def plot_tsat_vs_N(
    tsat_rows: list[tuple[str, float, int, float]],
    chains: list[str],
    lambdas: list[float],
    out_path: str = "fig_fss_tsat_vs_N.png",
) -> None:
    lams_nonzero = [lam for lam in lambdas if lam > 0.0]
    if not lams_nonzero:
        return
    cmap = plt.cm.plasma(np.linspace(0.05, 0.95, len(lams_nonzero)))
    style = {"fibonacci": "--", "tribonacci": "-"}
    marker = {"fibonacci": "^", "tribonacci": "o"}

    fig, ax = plt.subplots(figsize=(8, 5), dpi=140)
    for ch in chains:
        for i, lam in enumerate(lams_nonzero):
            pts = [(n, tsat) for (c, l, n, tsat) in tsat_rows
                   if c == ch and l == lam and np.isfinite(tsat)]
            if not pts:
                continue
            pts.sort()
            ns, tsats = zip(*pts)
            ax.loglog(ns, tsats, style.get(ch, "-") + marker.get(ch, "o"),
                      color=cmap[i], lw=1.4, ms=6,
                      label=f"{ch[:4]}  lambda={lam:.1f}")

    ax.set_xlabel("chain length  N")
    ax.set_ylabel("t_sat  (saturation crossover time)")
    ax.set_title("t_sat(N) -- finite-size saturation crossover")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(ncol=2, fontsize=8, loc="best")
    plt.tight_layout()
    plt.savefig(out_path, dpi=140)
    plt.close()
    print(f"  -> {out_path}")


def write_csvs(
    d2_fits: dict[str, tuple[float, float]],
    alpha_rows: list[tuple[str, float, int, float, int]],
    tsat_rows: list[tuple[str, float, int, float]],
) -> None:
    with open("fss_D2.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chain", "D2", "log_A"])
        for ch, (D2, log_A) in d2_fits.items():
            w.writerow([ch, f"{D2:.6f}", f"{log_A:.6f}"])
    print("  -> fss_D2.csv")

    with open("fss_alpha.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chain", "lambda", "N", "alpha", "n_points"])
        for ch, lam, n, alpha, n_pts in alpha_rows:
            w.writerow([ch, f"{lam:.4f}", n, f"{alpha:.6f}", n_pts])
    print("  -> fss_alpha.csv")

    with open("fss_tsat.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chain", "lambda", "N", "t_sat"])
        for ch, lam, n, tsat in tsat_rows:
            val = f"{tsat:.2f}" if np.isfinite(tsat) else "inf"
            w.writerow([ch, f"{lam:.4f}", n, val])
    print("  -> fss_tsat.csv")


# ---------------------------------------------------------------------------
# 5. Entry point
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Finite-size scaling analysis of long-time DNLS evolution."
    )
    ap.add_argument(
        "--n-values", type=int, nargs="+", default=[200, 500, 1000, 2000],
        help="chain lengths to include (default: 200 500 1000 2000)",
    )
    ap.add_argument(
        "--t-max", type=float, default=10000.0,
        help="maximum time to use from each CSV (default: 10000)",
    )
    ap.add_argument(
        "--lambda-sweep",
        default="ipr_lambda1p5_T1e4.csv",
        metavar="PATH",
        help=(
            "optional multi-N lambda-sweep CSV to merge into the analysis "
            "(produced by dnls_lambda1p5_sweep.py; default: ipr_lambda1p5_T1e4.csv). "
            "Skipped silently if the file does not exist."
        ),
    )
    ap.add_argument(
        "--no-plots", action="store_true",
        help="skip figure generation",
    )
    args = ap.parse_args()

    print(f"\nFinite-size scaling analysis  N={args.n_values}  T_max={args.t_max:.0f}")
    print("-" * 72)
    print("Loading data ...")
    all_data = load_all(args.n_values, args.t_max)

    # Merge optional lambda-sweep file (e.g. ipr_lambda1p5_T1e4.csv).
    # The file uses a multi-N format with an extra 'N' column.
    if os.path.exists(args.lambda_sweep):
        print(f"  merging lambda-sweep file: {args.lambda_sweep}")
        sweep_data = load_lambda_sweep(args.lambda_sweep, args.t_max)
        for n, nd in sweep_data.items():
            if n not in all_data:
                all_data[n] = {}
            all_data[n].update(nd)
        print(
            f"  merged N={sorted(sweep_data.keys())} "
            f"lambda={sorted({lam for ndd in sweep_data.values() for _, lam in ndd})}"
        )
    else:
        print(f"  [skip] lambda-sweep file not found: {args.lambda_sweep}")

    if not all_data:
        print("ERROR: no data loaded.", file=sys.stderr)
        return 1

    # infer chains and lambdas from loaded data
    chains = sorted({ch for nd in all_data.values() for ch, _ in nd.keys()})
    lambdas = sorted({lam for nd in all_data.values() for _, lam in nd.keys()})

    print(f"\nChains found : {chains}")
    print(f"Lambdas found: {lambdas}")
    print(f"N values loaded: {sorted(all_data.keys())}\n")

    d2_fits = analyze_D2(all_data, chains)
    alpha_rows = analyze_alpha(all_data, chains, lambdas)
    tsat_rows = analyze_tsat(all_data, chains, lambdas)

    print("Writing CSVs ...")
    write_csvs(d2_fits, alpha_rows, tsat_rows)

    if not args.no_plots:
        print("\nGenerating figures ...")
        plot_D2(all_data, chains, d2_fits)
        plot_alpha_vs_N(all_data, alpha_rows, chains, lambdas)
        plot_tsat_vs_N(tsat_rows, chains, lambdas)

    print("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
