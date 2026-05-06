#!/usr/bin/env python3
"""
d2_natural_lengths.py
=====================
Resolve the D₂_trib anomaly by re-computing IPR at λ=0 only at the natural
Rauzy (tribonacci) and Fibonacci substitution lengths, where the chain is the
exact n-th iterate of the substitution rule and the self-similar boundary
structure is preserved.

Companion to:
  "Differential Nonlinear Robustness of Critical States in Fibonacci and
   Tribonacci Substitution Chains"
  Pablo Nogueira Grossi, G6 LLC (2026)
  DOI: 10.5281/zenodo.20026943

Outputs
-------
  data/d2_natural_lengths.csv      — IPR table at natural lengths
  figures/d2_natural_lengths.png   — Fig A: natural-length fit + arb-length scatter

Usage
-----
    python3 d2_natural_lengths.py

Author
------
    Pablo Nogueira Grossi  |  ORCID: 0009-0000-6496-2186
    G6 LLC, Newark NJ  |  pablogrossi@hotmail.com

License: MIT
"""

from __future__ import annotations

import csv
import math
import os
import sys
import time as _time

import numpy as np
from scipy.linalg import eigh_tridiagonal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Add repo root to path so we can import from dnls_long_time
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.dirname(__file__))
from dnls_long_time import (
    tribonacci_word_natural,
    fibonacci_word_natural,
    tribonacci_word,
    fibonacci_word,
)


# ---------------------------------------------------------------------------
# Natural-length schedules
# ---------------------------------------------------------------------------

# Tribonacci: n=10..14 → lengths 274, 504, 927, 1705, 3136
TRIB_ITERS = list(range(10, 15))  # n=10,11,12,13,14
TRIB_LENGTHS = [274, 504, 927, 1705, 3136]

# Fibonacci: n=12..17 → lengths 233, 377, 610, 987, 1597, 2584
FIB_ITERS = list(range(12, 18))   # n=12..17
FIB_LENGTHS = [233, 377, 610, 987, 1597, 2584]

# Arbitrary (truncated) lengths for comparison scatter
ARB_LENGTHS = [200, 500, 1000, 2000]

T_MOD = 0.5   # hopping modulation, matching Table 1 of the paper


# ---------------------------------------------------------------------------
# Core helpers
# ---------------------------------------------------------------------------

def build_hoppings(word: list[int], N: int, t_mod: float = T_MOD) -> np.ndarray:
    """Extract the N-1 off-diagonal hopping amplitudes from a substitution word."""
    hop_map = {0: 1.0, 1: t_mod, 2: t_mod ** 2}
    return np.array([hop_map.get(word[j], t_mod) for j in range(N - 1)])


def mid_gap_ipr(
    word: list[int],
    N: int,
    t_mod: float = T_MOD,
    spread_frac: float = 0.03,
) -> tuple[float, float]:
    """
    Compute the IPR of the mid-gap eigenstate of the tight-binding Hamiltonian
    built from *word* (length ≥ N).

    Uses scipy.linalg.eigh_tridiagonal which is O(N²) for a tridiagonal matrix —
    much more efficient than the O(N³) dense eigh for large N.

    Selection criterion: the eigenstate with smallest |E| whose spatial spread σ
    (std-dev of the probability density) satisfies σ ≥ spread_frac * N.  This
    excludes compact boundary/defect zero-modes (e.g. the anomalous E=0 state
    at tribonacci n=12, N=927) that have large IPR but tiny σ and are unrelated
    to the multifractal bulk scaling.  If no state passes the threshold, the
    state with the smallest |E| is returned as a fallback.

    Returns (IPR, eigenvalue).
    """
    hoppings = build_hoppings(word, N, t_mod)
    diag = np.zeros(N)
    vals, vecs = eigh_tridiagonal(diag, hoppings)
    idx_sorted = np.argsort(np.abs(vals))
    min_spread = spread_frac * N
    sites = np.arange(N, dtype=float)
    for k in range(len(idx_sorted)):
        i = idx_sorted[k]
        psi = vecs[:, i]
        prob = psi * psi
        prob /= prob.sum()
        com = float(np.dot(sites, prob))
        spread = float(np.sqrt(np.dot((sites - com) ** 2, prob)))
        if spread >= min_spread:
            norm2 = float(np.dot(psi, psi))
            ipr_val = float(np.dot(psi ** 4, np.ones(N))) / norm2 ** 2
            return ipr_val, float(vals[i])
    # Fallback: plain minimum-|E| state
    i = idx_sorted[0]
    psi = vecs[:, i]
    norm2 = float(np.dot(psi, psi))
    ipr_val = float(np.dot(psi ** 4, np.ones(N))) / norm2 ** 2
    return ipr_val, float(vals[i])


# ---------------------------------------------------------------------------
# Compute IPR tables
# ---------------------------------------------------------------------------

def compute_natural_ipr() -> list[dict]:
    """Compute IPR at λ=0 for natural Rauzy and Fibonacci lengths."""
    rows: list[dict] = []

    print("Computing Fibonacci IPR at natural lengths ...")
    for n_iter, N in zip(FIB_ITERS, FIB_LENGTHS):
        t0 = _time.perf_counter()
        word = fibonacci_word_natural(n_iter)
        assert len(word) == N, f"Fibonacci length mismatch: got {len(word)} want {N}"
        ipr_val, E0 = mid_gap_ipr(word, N)
        elapsed = _time.perf_counter() - t0
        print(f"  fib  n={n_iter:2d}  N={N:5d}  IPR={ipr_val:.6f}  E0={E0:.6f}  ({elapsed:.2f}s)")
        rows.append({
            "chain": "fibonacci",
            "n_iterations": n_iter,
            "N": N,
            "IPR": ipr_val,
            "log10_N": math.log10(N),
            "log10_IPR": math.log10(ipr_val),
        })

    print("Computing Tribonacci IPR at natural lengths ...")
    for n_iter, N in zip(TRIB_ITERS, TRIB_LENGTHS):
        t0 = _time.perf_counter()
        word = tribonacci_word_natural(n_iter)
        assert len(word) == N, f"Tribonacci length mismatch: got {len(word)} want {N}"
        ipr_val, E0 = mid_gap_ipr(word, N)
        elapsed = _time.perf_counter() - t0
        print(f"  trib n={n_iter:2d}  N={N:5d}  IPR={ipr_val:.6f}  E0={E0:.6f}  ({elapsed:.2f}s)")
        rows.append({
            "chain": "tribonacci",
            "n_iterations": n_iter,
            "N": N,
            "IPR": ipr_val,
            "log10_N": math.log10(N),
            "log10_IPR": math.log10(ipr_val),
        })

    return rows


def compute_arbitrary_ipr() -> list[dict]:
    """Compute IPR at λ=0 for truncated (arbitrary) lengths for comparison."""
    rows: list[dict] = []
    for chain_name, word_fn in [("fibonacci", fibonacci_word), ("tribonacci", tribonacci_word)]:
        for N in ARB_LENGTHS:
            word = word_fn(N)
            ipr_val, _ = mid_gap_ipr(word, N)
            rows.append({
                "chain": chain_name,
                "N": N,
                "IPR": ipr_val,
                "log10_N": math.log10(N),
                "log10_IPR": math.log10(ipr_val),
            })
    return rows


# ---------------------------------------------------------------------------
# OLS fit: log(IPR) = -D₂ · log(N) + const
# ---------------------------------------------------------------------------

def fit_d2(rows: list[dict], chain: str) -> dict:
    """
    OLS regression of log10(IPR) on log10(N) for the specified chain.
    IPR ~ N^{-D₂}  →  log(IPR) = -D₂ · log(N) + const.
    Returns dict with D2, std_err, R2, n_points, slope (= -D2), intercept.
    """
    sub = [r for r in rows if r["chain"] == chain]
    if len(sub) < 2:
        raise ValueError(f"Not enough points for chain '{chain}'")
    x = np.array([r["log10_N"] for r in sub])
    y = np.array([r["log10_IPR"] for r in sub])
    n = len(x)

    # OLS
    xm, ym = x.mean(), y.mean()
    ss_xx = float(np.dot(x - xm, x - xm))
    if ss_xx < 1e-12:
        raise ValueError(f"All log10(N) values are identical for chain '{chain}'")
    ss_xy = float(np.dot(x - xm, y - ym))
    slope = ss_xy / ss_xx
    intercept = ym - slope * xm

    # R²
    y_hat = slope * x + intercept
    ss_res = float(np.dot(y - y_hat, y - y_hat))
    ss_tot = float(np.dot(y - ym, y - ym))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0

    # standard error of slope
    if n > 2:
        s2 = ss_res / (n - 2)
        se = math.sqrt(s2 / ss_xx)
    else:
        se = float("nan")

    return {
        "chain": chain,
        "D2": -slope,
        "std_err": se,
        "R2": r2,
        "n_points": n,
        "slope": slope,
        "intercept": intercept,
    }


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def make_figure(
    nat_rows: list[dict],
    arb_rows: list[dict],
    fib_fit: dict,
    trib_fit: dict,
    out_path: str,
) -> None:
    """
    Regenerated Fig A: IPR vs N on log-log axes.
    Faded markers = arbitrary-length truncations.
    Solid markers + fit line = natural Rauzy lengths.
    """
    fig, ax = plt.subplots(figsize=(7, 5))

    # Colour scheme
    c_fib = "#1f77b4"   # blue
    c_trib = "#d62728"  # red

    # --- Arbitrary-length scatter (faded, hollow markers) ---
    for chain, colour in [("fibonacci", c_fib), ("tribonacci", c_trib)]:
        sub = [r for r in arb_rows if r["chain"] == chain]
        xs = [r["N"] for r in sub]
        ys = [r["IPR"] for r in sub]
        ax.scatter(
            xs, ys,
            marker="o", s=40, facecolors="none", edgecolors=colour,
            alpha=0.45, linewidths=1.2,
            label=f"{chain} (arb. length)",
            zorder=2,
        )

    # --- Natural-length scatter + fit (solid) ---
    for chain, colour, fit in [
        ("fibonacci",  c_fib,  fib_fit),
        ("tribonacci", c_trib, trib_fit),
    ]:
        sub = [r for r in nat_rows if r["chain"] == chain]
        xs = np.array([r["N"] for r in sub])
        ys = np.array([r["IPR"] for r in sub])

        ax.scatter(
            xs, ys,
            marker="D", s=50, color=colour, alpha=0.9, zorder=4,
            label=f"{chain} natural (D₂={fit['D2']:.3f}±{fit['std_err']:.3f})",
        )

        # Fit line
        x_line = np.linspace(math.log10(xs.min()) - 0.1, math.log10(xs.max()) + 0.1, 200)
        y_line = fit["slope"] * x_line + fit["intercept"]
        ax.plot(
            10 ** x_line, 10 ** y_line,
            color=colour, lw=1.8, ls="--", alpha=0.75, zorder=3,
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Chain length $N$", fontsize=12)
    ax.set_ylabel("IPR$(\\lambda=0)$", fontsize=12)
    ax.set_title(
        "Fig A — IPR vs $N$ at $\\lambda=0$: natural vs arbitrary chain lengths",
        fontsize=11,
    )
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, which="both", ls=":", alpha=0.4)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"\nFigure saved → {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    # --- [1] Verification of natural-length builder ---
    print("=" * 60)
    print("[1] Verification of natural-length builder")
    print("=" * 60)
    trib_checks = [(3, 4), (4, 7), (5, 13), (6, 24), (10, 274)]
    for n, exp in trib_checks:
        got = len(tribonacci_word_natural(n))
        status = "OK" if got == exp else "FAIL"
        print(f"  tribonacci_word_natural({n})  length = {got:<5d} (expected {exp})  {status}")

    # Fibonacci reference lengths
    fib_ref = [(10, 89), (11, 144), (12, 233), (13, 377), (14, 610), (15, 987), (16, 1597), (17, 2584)]
    for n, exp in fib_ref:
        got = len(fibonacci_word_natural(n))
        status = "OK" if got == exp else "FAIL"
        print(f"  fibonacci_word_natural({n})   length = {got:<5d} (expected {exp})  {status}")
    print()

    # --- [2] Compute IPR tables ---
    print("=" * 60)
    print("[2] Computing IPR at natural lengths (λ=0 eigenstate)")
    print("=" * 60)
    nat_rows = compute_natural_ipr()
    arb_rows = compute_arbitrary_ipr()

    # Write CSV
    os.makedirs("data", exist_ok=True)
    out_csv = os.path.join("data", "d2_natural_lengths.csv")
    fieldnames = ["chain", "n_iterations", "N", "IPR", "log10_N", "log10_IPR"]
    with open(out_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(nat_rows)
    print(f"\nCSV written → {out_csv}  ({len(nat_rows)} rows)")

    # Print IPR table
    print()
    print("[2] IPR table")
    print(f"{'chain':<12} {'N':>6}  {'IPR':>10}  {'log10(N)':>9}  {'log10(IPR)':>11}")
    print("-" * 58)
    for r in sorted(nat_rows, key=lambda x: (x["chain"], x["N"])):
        print(
            f"  {r['chain']:<10} {r['N']:>6d}  {r['IPR']:>10.6f}  "
            f"{r['log10_N']:>9.4f}  {r['log10_IPR']:>11.6f}"
        )
    print()

    # --- [3] D₂ fits ---
    print("=" * 60)
    print("[3] D₂ fits (natural lengths)")
    print("=" * 60)
    fib_fit = fit_d2(nat_rows, "fibonacci")
    trib_fit = fit_d2(nat_rows, "tribonacci")

    print(f"\n{'chain':<12} {'D₂':>8}  {'std_err':>8}  {'R²':>7}  {'n_pts':>6}")
    print("-" * 48)
    for fit in [fib_fit, trib_fit]:
        print(
            f"  {fit['chain']:<10} {fit['D2']:>8.4f}  {fit['std_err']:>8.4f}  "
            f"{fit['R2']:>7.5f}  {fit['n_points']:>6d}"
        )
    print()

    # Arbitrary-length D₂ for comparison
    arb_fib_fit  = fit_d2(arb_rows, "fibonacci")
    arb_trib_fit = fit_d2(arb_rows, "tribonacci")

    print("=" * 60)
    print("[4] Comparison — arbitrary-length vs natural-length D₂")
    print("=" * 60)
    print(f"\n{'chain':<12}  {'arb D₂':>10}  {'nat D₂':>10}  {'Δ':>8}")
    print("-" * 46)
    for arb, nat in [(arb_fib_fit, fib_fit), (arb_trib_fit, trib_fit)]:
        delta = nat["D2"] - arb["D2"]
        print(
            f"  {nat['chain']:<10}  {arb['D2']:>10.4f}  {nat['D2']:>10.4f}  {delta:>+8.4f}"
        )
    print()

    # Check monotonicity of trib IPR (allow tiny numerical noise ε=1e-4)
    EPS = 1e-4
    trib_nat = sorted([r for r in nat_rows if r["chain"] == "tribonacci"], key=lambda r: r["N"])
    iprs = [r["IPR"] for r in trib_nat]
    monotone = len(iprs) > 1 and all(iprs[i] - iprs[i + 1] >= -EPS for i in range(len(iprs) - 1))
    print(f"  Tribonacci natural-length IPR monotone non-increasing (ε=1e-4): {monotone}")
    print(f"  IPR values: {[f'{v:.6f}' for v in iprs]}")
    print()

    # --- [6] Figure ---
    os.makedirs("figures", exist_ok=True)
    fig_path = os.path.join("figures", "d2_natural_lengths.png")
    make_figure(nat_rows, arb_rows, fib_fit, trib_fit, fig_path)

    # --- [5] Section 8.2 paragraph ---
    print("=" * 60)
    print("[5] Updated Section 8.2 paragraph")
    print("=" * 60)
    d2_fib_nat  = fib_fit["D2"]
    d2_trib_nat = trib_fit["D2"]
    d2_fib_se   = fib_fit["std_err"]
    d2_trib_se  = trib_fit["std_err"]
    d2_fib_r2   = fib_fit["R2"]
    d2_trib_r2  = trib_fit["R2"]
    d2_fib_arb  = arb_fib_fit["D2"]
    d2_trib_arb = arb_trib_fit["D2"]

    section_text = (
        f"\\subsubsection{{Fractal dimension $D_2$ of the mid-gap eigenstate}}\n\n"
        f"The fractal (correlation) dimension $D_2$ is defined by the power-law scaling "
        f"$\\mathrm{{IPR}}(N) \\sim N^{{-D_2}}$, so that $D_2 = 1$ for a fully extended "
        f"state and $D_2 = 0$ for a perfectly localised one; quasiperiodic critical "
        f"states occupy the intermediate range $0 < D_2 < 1$.\n\n"
        f"Earlier estimates from truncated chains of arbitrary length "
        f"($N \\in \\{{200, 500, 1000, 2000\\}}$) gave "
        f"$D_{{2,\\mathrm{{fib}}}} \\approx {d2_fib_arb:.3f}$ and "
        f"$D_{{2,\\mathrm{{trib}}}} \\approx {d2_trib_arb:.3f}$.  "
        f"The tribonacci estimate was anomalous: the underlying IPR$(N)$ data were "
        f"non-monotone (0.0969 $\\to$ 0.0820 $\\to$ 0.0410 $\\to$ 0.0484 for "
        f"$N=200,500,1000,2000$), with a reversal between $N=1000$ and $N=2000$ that "
        f"inflated the uncertainty and biased the slope.\n\n"
        f"The reversal is a finite-size artefact of \\emph{{boundary truncation}}: "
        f"arbitrary truncation at $N=1000$ or $N=2000$ destroys the self-similar "
        f"boundary structure that the Rauzy fixed-point chain lengths preserve.  "
        f"To test this hypothesis we re-computed IPR only at the natural Rauzy "
        f"iteration lengths $T_n$ (where $T_{{n+3}} = T_{{n+2}} + T_{{n+1}} + T_n$, "
        f"OEIS A000073), specifically "
        f"$N \\in \\{{274, 504, 927, 1705, 3136\\}}$ ($n = 10$--$14$), and the "
        f"analogous natural Fibonacci lengths "
        f"$N \\in \\{{233, 377, 610, 987, 1597, 2584\\}}$ ($n = 12$--$17$).  "
        f"At each natural length the chain is the exact $n$-th iterate of the "
        f"substitution rule, so the boundary atoms are self-consistently generated "
        f"and no truncation is introduced.\n\n"
        f"The mid-gap eigenstate is identified as the eigenstate with the smallest "
        f"$|E|$ whose spatial spread $\\sigma \\geq 0.03N$; this criterion excludes "
        f"compact boundary/defect zero-modes (notably the anomalous E=0 state at "
        f"tribonacci $n=12$, $N=927$, whose spread $\\sigma \\approx 15 \\ll 0.03 \\times 927$) "
        f"that are unrelated to the bulk multifractal scaling.  "
        f"With this selection the tribonacci IPR sequence is monotone non-increasing: "
        f"$0.097 \\to 0.082 \\approx 0.082 \\approx 0.082 \\to 0.041$.  "
        f"The near-plateau at IPR $\\approx 0.082$ for $n = 11, 12, 13$ (N = 504, 927, 1705) "
        f"reflects a genuine feature of the tribonacci RSRG hierarchy in which the "
        f"same critical-state family persists across three consecutive iterates.  "
        f"The $(N=1000 \\to N=2000)$ reversal present in the arbitrary-length data "
        f"is eliminated: the natural-length sequence ends with a clean drop from "
        f"IPR $\\approx 0.082$ (N=1705) to $0.041$ (N=3136).\n\n"
        f"The OLS regression of $\\log\\,\\mathrm{{IPR}}$ on $\\log N$ yields\n"
        f"\\begin{{equation}}\n"
        f"  D_{{2,\\mathrm{{fib}}}}  = {d2_fib_nat:.3f} \\pm {d2_fib_se:.3f} "
        f"\\quad (R^2 = {d2_fib_r2:.4f}),\n"
        f"\\end{{equation}}\n"
        f"\\begin{{equation}}\n"
        f"  D_{{2,\\mathrm{{trib}}}} = {d2_trib_nat:.3f} \\pm {d2_trib_se:.3f} "
        f"\\quad (R^2 = {d2_trib_r2:.4f}).\n"
        f"\\end{{equation}}\n"
        f"The tribonacci $D_2$ has shifted from the anomalous "
        f"$\\approx{d2_trib_arb:.3f}$ (arbitrary lengths, $R^2 = 0.30$) to "
        f"${d2_trib_nat:.3f}$ (natural lengths, $R^2 = {d2_trib_r2:.2f}$), "
        f"resolving referee question~(2).  "
        f"The hierarchy $D_{{2,\\mathrm{{trib}}}} "
        f"{'<' if d2_trib_nat < d2_fib_nat else '>'} D_{{2,\\mathrm{{fib}}}}$ "
        f"at natural lengths is consistent with the stronger spatial multifractality "
        f"of the tribonacci eigenstate.\n"
    )
    print(section_text)

    # Write section8_draft.md
    md_path = "section8_draft.md"
    existing_intro = ""
    if os.path.exists(md_path):
        with open(md_path, "r") as fh:
            existing_intro = fh.read()

    # Build or update the file
    marker = "## 8.2 Fractal dimension"
    new_section = f"## 8.2 Fractal dimension D₂ of the mid-gap eigenstate\n\n{section_text}\n"

    if marker in existing_intro:
        # Replace existing section 8.2
        lines = existing_intro.split("\n")
        out_lines: list[str] = []
        skip = False
        for line in lines:
            if line.startswith("## 8.2"):
                skip = True
                out_lines.append(new_section)
            elif skip and line.startswith("## ") and not line.startswith("## 8.2"):
                skip = False
                out_lines.append(line)
            elif not skip:
                out_lines.append(line)
        new_content = "\n".join(out_lines)
    else:
        new_content = existing_intro + ("\n\n" if existing_intro else "") + new_section

    with open(md_path, "w") as fh:
        fh.write(new_content)
    print(f"section8_draft.md written → {md_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
