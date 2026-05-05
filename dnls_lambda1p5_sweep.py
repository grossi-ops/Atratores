#!/usr/bin/env python3
"""
dnls_lambda1p5_sweep.py
=======================
Fine-grained nonlinearity sweep around the λ ≈ 1.5 crossover on Fibonacci
and Tribonacci substitution chains.

Context
-------
`dnls_nbonacci.py` (Table 1 of the companion paper) identifies λ = 1.5 as a
key transition point where the IPR of the Fibonacci mid-gap state has already
dropped substantially while the Tribonacci state remains more localised.  This
script resolves that crossover region with step Δλ = 0.1 (configurable) and
records both the short-time IPR (after a fixed evolution time T) and the
instantaneous linear-limit IPR for reference.

Companion to
------------
  "Differential Nonlinear Robustness of Critical States in Fibonacci and
   Tribonacci Substitution Chains"
  Pablo Nogueira Grossi, G6 LLC (2026)
  DOI (this paper): 10.5281/zenodo.20026943
  DOI (Vol. I):     10.5281/zenodo.19117400

Outputs
-------
  lambda1p5_sweep.csv       — IPR vs λ, long-format (lambda, chain, IPR, norm)
  lambda1p5_sweep_fig.png   — IPR(λ) plot for both chains with crossover detail

Usage
-----
    python3 dnls_lambda1p5_sweep.py                   # default parameters
    python3 dnls_lambda1p5_sweep.py --lam-min 0.5 --lam-max 2.5 --lam-step 0.05
    python3 dnls_lambda1p5_sweep.py --no-plot --quiet

Dependencies
------------
    numpy, scipy, matplotlib  (standard scientific Python)

Author
------
    Pablo Nogueira Grossi  |  ORCID: 0009-0000-6496-2186
    G6 LLC, Newark NJ  |  pablogrossi@hotmail.com
    GitHub: https://github.com/TOTOGT/AXLE

License: MIT
"""

from __future__ import annotations

import argparse
import csv
import sys
import time as _time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigh

# ---------------------------------------------------------------------------
# Default sweep parameters
# ---------------------------------------------------------------------------

N_SITES      = 500      # chain length (matches Table 1 of the paper)
T_MOD        = 0.5      # hopping modulation (generic incommensurate value)
T_EVO        = 50.0     # integration time (matches dnls_nbonacci.py Table 1)
LAM_MIN      = 0.0      # sweep start
LAM_MAX      = 3.0      # sweep end  (well past the λ=1.5 crossover)
LAM_STEP     = 0.1      # step size  (10× finer than dnls_nbonacci.py)
OUT_CSV      = "lambda1p5_sweep.csv"
OUT_FIG      = "lambda1p5_sweep_fig.png"
RTOL         = 1e-7
ATOL         = 1e-9
NORM_TOL     = 1e-4     # warn if |L2-norm| drifts more than this


# ---------------------------------------------------------------------------
# 1.  Substitution words (verbatim from dnls_nbonacci.py)
# ---------------------------------------------------------------------------

def fibonacci_word(length: int) -> list[int]:
    """Fibonacci substitution: A→AB, B→A (encoded 0,1)."""
    word = [0]
    rules = {0: [0, 1], 1: [0]}
    while len(word) < length:
        word = [s for c in word for s in rules[c]]
    return word[:length]


def tribonacci_word(length: int) -> list[int]:
    """Rauzy (tribonacci) substitution: A→AB, B→AC, C→A (encoded 0,1,2)."""
    word = [0]
    rules = {0: [0, 1], 1: [0, 2], 2: [0]}
    while len(word) < length:
        word = [s for c in word for s in rules[c]]
    return word[:length]


# ---------------------------------------------------------------------------
# 2.  Tight-binding Hamiltonian and linear-limit utilities
# ---------------------------------------------------------------------------

def build_hamiltonian(
    word: list[int], N: int, t_mod: float = 0.5
) -> tuple[np.ndarray, np.ndarray]:
    """Tridiagonal tight-binding Hamiltonian on N sites."""
    hop_map = {0: 1.0, 1: t_mod, 2: t_mod ** 2}
    hoppings = np.array([hop_map.get(word[j], t_mod) for j in range(N - 1)])
    H = np.zeros((N, N))
    for j in range(N - 1):
        H[j, j + 1] = hoppings[j]
        H[j + 1, j] = hoppings[j]
    return H, hoppings


def mid_gap_state(H: np.ndarray) -> tuple[np.ndarray, float]:
    """Eigenstate of H whose eigenvalue is closest to E = 0."""
    vals, vecs = eigh(H)
    idx = np.argmin(np.abs(vals))
    return vecs[:, idx], float(vals[idx])


def ipr(psi: np.ndarray) -> float:
    """Inverse participation ratio: Σ|ψ_j|⁴ / (Σ|ψ_j|²)²."""
    norm2 = float(np.sum(np.abs(psi) ** 2))
    return float(np.sum(np.abs(psi) ** 4)) / norm2 ** 2


# ---------------------------------------------------------------------------
# 3.  DNLS right-hand side and integrator (matching dnls_long_time.py)
# ---------------------------------------------------------------------------

def dnls_rhs(
    t: float,
    state: np.ndarray,
    lam: float,
    hoppings: np.ndarray,
) -> np.ndarray:
    """
    Real-valued formulation of DNLS:
        i dψ_j/dt = −Σ_{j'} H_{jj'} ψ_{j'} + λ|ψ_j|² ψ_j

    state = [Re(ψ); Im(ψ)], length 2N.
    """
    N = len(state) >> 1
    x = state[:N]
    y = state[N:]

    Hx = np.zeros(N)
    Hx[:-1] += hoppings * x[1:]
    Hx[1:]  += hoppings * x[:-1]

    Hy = np.zeros(N)
    Hy[:-1] += hoppings * y[1:]
    Hy[1:]  += hoppings * y[:-1]

    nl = x * x + y * y          # |ψ_j|²
    dxdt = -Hy + lam * nl * y
    dydt =  Hx - lam * nl * x
    return np.concatenate([dxdt, dydt])


def evolve_dnls(
    psi0: np.ndarray,
    hoppings: np.ndarray,
    lam: float,
    T: float = 50.0,
    rtol: float = RTOL,
    atol: float = ATOL,
) -> tuple[np.ndarray, float]:
    """
    Integrate DNLS from t=0 to t=T with RK45 and return the final state.

    Returns
    -------
    psi_final : (N,) complex ndarray
    norm_final : float  (should remain ≈ initial norm)
    """
    N = len(psi0)
    # Normalise to unit L2 norm
    psi0 = psi0 / np.sqrt(np.dot(psi0, psi0))
    state0 = np.concatenate([psi0, np.zeros_like(psi0)])

    sol = solve_ivp(
        dnls_rhs,
        [0.0, T],
        state0,
        args=(lam, hoppings),
        method="RK45",
        rtol=rtol,
        atol=atol,
        max_step=0.1,
        dense_output=False,
    )
    zf = sol.y[:, -1]
    psi_f = zf[:N] + 1j * zf[N:]
    return psi_f, float(np.sqrt(np.sum(np.abs(psi_f) ** 2)))


# ---------------------------------------------------------------------------
# 4.  Lambda sweep
# ---------------------------------------------------------------------------

def lambda_sweep(
    n: int = N_SITES,
    t_mod: float = T_MOD,
    t_evo: float = T_EVO,
    lam_min: float = LAM_MIN,
    lam_max: float = LAM_MAX,
    lam_step: float = LAM_STEP,
    out_csv: str = OUT_CSV,
    out_fig: str = OUT_FIG,
    no_plot: bool = False,
    verbose: bool = True,
) -> list[dict]:
    """
    Sweep λ from lam_min to lam_max in steps of lam_step.

    For each λ and each chain (Fibonacci, Tribonacci):
      - evolve the mid-gap eigenstate under DNLS for time t_evo
      - record IPR of the final state and the L2-norm (conservation check)

    Outputs
    -------
    CSV with columns: lambda, chain, IPR_linear, IPR_final, norm_final
    PNG figure of IPR_final vs λ for both chains (unless no_plot=True)

    Returns the list of result dicts.
    """
    # Build λ array — round to avoid floating-point tick accumulation
    n_steps = round((lam_max - lam_min) / lam_step)
    lambdas = [round(lam_min + i * lam_step, 10) for i in range(n_steps + 1)]

    # ── Build chains once ────────────────────────────────────────────────────
    if verbose:
        print(f"Building chains  N={n}  t_mod={t_mod} ...")

    # n + 1 so that the Hamiltonian builder can read n bond letters (indices 0..n-1)
    # from a word of at least length n + 1.
    chains = {
        "fibonacci":  fibonacci_word(n + 1),
        "tribonacci": tribonacci_word(n + 1),
    }
    hamiltonians: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    psi0s: dict[str, tuple[np.ndarray, float]] = {}
    ipr0s: dict[str, float] = {}

    for name, word in chains.items():
        H, hops = build_hamiltonian(word, n, t_mod)
        psi0, E0 = mid_gap_state(H)
        hamiltonians[name] = (H, hops)
        psi0s[name] = (psi0, E0)
        ipr0s[name] = ipr(psi0)
        if verbose:
            print(f"  {name:>12}  E0={E0:.6f}  IPR(linear)={ipr0s[name]:.6f}")

    if verbose:
        print()

    # ── Sweep ────────────────────────────────────────────────────────────────
    n_runs = len(lambdas) * len(chains)
    run_idx = 0
    rows: list[dict] = []

    header = (
        f"{'lambda':>8} | {'chain':>12} | {'IPR_linear':>12} | "
        f"{'IPR_final':>12} | {'norm_final':>12} | {'elapsed_s':>9}"
    )
    if verbose:
        print(header)
        print("-" * len(header))

    for lam in lambdas:
        for name in ("fibonacci", "tribonacci"):
            run_idx += 1
            _, hops = hamiltonians[name]
            psi0, _  = psi0s[name]

            t0 = _time.perf_counter()
            psi_f, norm_f = evolve_dnls(psi0, hops, lam, T=t_evo)
            elapsed = _time.perf_counter() - t0

            ipr_f = ipr(psi_f)
            norm_warn = "" if abs(norm_f - 1.0) <= NORM_TOL else "  *** NORM LEAK ***"

            row = {
                "lambda":     lam,
                "chain":      name,
                "IPR_linear": ipr0s[name],
                "IPR_final":  ipr_f,
                "norm_final": norm_f,
            }
            rows.append(row)

            if verbose:
                print(
                    f"{lam:>8.2f} | {name:>12} | {ipr0s[name]:>12.6f} | "
                    f"{ipr_f:>12.6f} | {norm_f:>12.6f} | {elapsed:>9.1f}s"
                    f"{norm_warn}"
                )

    # ── Write CSV ────────────────────────────────────────────────────────────
    fields = ["lambda", "chain", "IPR_linear", "IPR_final", "norm_final"]
    with open(out_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    if verbose:
        print(f"\nWrote {len(rows)} rows -> {out_csv}")

    # ── Figure ───────────────────────────────────────────────────────────────
    if not no_plot:
        _plot_sweep(rows, lambdas, out_fig, t_evo, n, verbose)

    return rows


# ---------------------------------------------------------------------------
# 5.  Plot helper
# ---------------------------------------------------------------------------

def _plot_sweep(
    rows: list[dict],
    lambdas: list[float],
    out_fig: str,
    t_evo: float,
    n: int,
    verbose: bool,
) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        if verbose:
            print("matplotlib not available — skipping plot")
        return

    # Separate chains
    fib_ipr  = [r["IPR_final"]  for r in rows if r["chain"] == "fibonacci"]
    trib_ipr = [r["IPR_final"]  for r in rows if r["chain"] == "tribonacci"]
    fib_lam  = [r["lambda"]     for r in rows if r["chain"] == "fibonacci"]
    trib_lam = [r["lambda"]     for r in rows if r["chain"] == "tribonacci"]
    fib_lin  = [r["IPR_linear"] for r in rows if r["chain"] == "fibonacci"]
    trib_lin = [r["IPR_linear"] for r in rows if r["chain"] == "tribonacci"]

    COL_FIB  = "#2166ac"
    COL_TRIB = "#d6604d"
    COL_GOLD = "#c9a84c"

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 7), sharex=True,
                                   gridspec_kw={"hspace": 0.05})

    # ── Top: IPR(λ) ──────────────────────────────────────────────────────────
    ax1.plot(fib_lam,  fib_ipr,  "s-",  color=COL_FIB,  lw=1.6, ms=4,
             label=f"Fibonacci (N={n})")
    ax1.plot(trib_lam, trib_ipr, "o-",  color=COL_TRIB, lw=1.6, ms=4,
             label=f"Tribonacci (N={n})")
    ax1.axhline(fib_lin[0],  color=COL_FIB,  lw=0.7, ls=":", alpha=0.6,
                label="Linear IPR (Fibonacci)")
    ax1.axhline(trib_lin[0], color=COL_TRIB, lw=0.7, ls=":", alpha=0.6,
                label="Linear IPR (Tribonacci)")
    ax1.axvline(1.5, color=COL_GOLD, lw=1.0, ls="--", alpha=0.8,
                label="λ = 1.5")
    ax1.set_ylabel("IPR after T = {:.0f}".format(t_evo))
    ax1.set_title(
        "Fine λ-sweep around the λ ≈ 1.5 crossover\n"
        "(Fibonacci vs Tribonacci substitution chains)"
    )
    ax1.legend(fontsize=8, loc="upper right", ncol=2)
    ax1.grid(True, alpha=0.3)

    # ── Bottom: IPR ratio tribonacci / fibonacci ──────────────────────────────
    lam_arr  = np.array(fib_lam)
    fib_arr  = np.array(fib_ipr)
    trib_arr = np.array(trib_ipr)
    ratio    = trib_arr / fib_arr

    ax2.plot(lam_arr, ratio, "D-", color=COL_GOLD, lw=1.6, ms=4,
             label="IPR ratio trib / fib")
    ax2.axhline(1.0, color="#555555", lw=0.8, ls=":")
    ax2.axvline(1.5, color=COL_GOLD, lw=1.0, ls="--", alpha=0.8)
    ax2.fill_between(lam_arr, ratio, 1.0, where=(ratio > 1),
                     alpha=0.12, color=COL_GOLD)
    ax2.set_xlabel("Nonlinearity strength λ")
    ax2.set_ylabel("IPR ratio (trib / fib)")
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_fig, dpi=150)
    plt.close()
    if verbose:
        print(f"Wrote figure -> {out_fig}")


# ---------------------------------------------------------------------------
# 6.  CLI entry point
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Fine-grained λ sweep around the λ≈1.5 crossover on "
            "Fibonacci and Tribonacci substitution chains."
        )
    )
    ap.add_argument("-N", "--sites",    type=int,   default=N_SITES,
                    help=f"chain length (default: {N_SITES})")
    ap.add_argument("--t-mod",          type=float, default=T_MOD,
                    help=f"hopping modulation t_mod (default: {T_MOD})")
    ap.add_argument("-T", "--t-evo",    type=float, default=T_EVO,
                    help=f"DNLS evolution time (default: {T_EVO})")
    ap.add_argument("--lam-min",        type=float, default=LAM_MIN,
                    help=f"sweep start λ (default: {LAM_MIN})")
    ap.add_argument("--lam-max",        type=float, default=LAM_MAX,
                    help=f"sweep end λ (default: {LAM_MAX})")
    ap.add_argument("--lam-step",       type=float, default=LAM_STEP,
                    help=f"sweep step Δλ (default: {LAM_STEP})")
    ap.add_argument("--out",            default=OUT_CSV,
                    help=f"output CSV path (default: {OUT_CSV})")
    ap.add_argument("--fig",            default=OUT_FIG,
                    help=f"output figure path (default: {OUT_FIG})")
    ap.add_argument("--no-plot",        action="store_true",
                    help="skip matplotlib figure")
    ap.add_argument("--quiet",          action="store_true",
                    help="suppress progress output")
    args = ap.parse_args()

    lambda_sweep(
        n        = args.sites,
        t_mod    = args.t_mod,
        t_evo    = args.t_evo,
        lam_min  = args.lam_min,
        lam_max  = args.lam_max,
        lam_step = args.lam_step,
        out_csv  = args.out,
        out_fig  = args.fig,
        no_plot  = args.no_plot,
        verbose  = not args.quiet,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
