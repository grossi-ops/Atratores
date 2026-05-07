#!/usr/bin/env python3
"""
dnls_tetrabonacci.py
====================
Tetrabonacci substitution chain: T=50, T=10^4 FSS, and D₂.

Companion to the follow-up paper extending the fibonacci → tribonacci
differential-nonlinear-robustness analysis to the tetrabonacci case.

Substitution rule: A -> AB, B -> AC, C -> AD, D -> A.
Hopping convention (geometric, ratio 0.5): t_A=1.0, t_B=0.5, t_C=0.25, t_D=0.125.
This convention is a choice; it is documented explicitly so referees can
distinguish it from alternative parameterisations.

Perron–Frobenius eigenvalue (tetrabonacci constant): largest root of
    x^4 - x^3 - x^2 - x - 1 = 0  ≈ 1.92756.

Tasks performed
---------------
B  Verify natural iteration lengths against OEIS A000078.
C  T=50, N=500, RK45, RTOL=1e-6 across all LAMBDAS  ->  data/tetrabonacci_T50_N500.csv
D  T=10^4, lambda=1.5, N in {500,1000,2000}, DOP853, RTOL=1e-9  ->  data/tetrabonacci_lambda1p5_T1e4.csv
E  D2 from mid-gap eigenstate IPR at natural lengths  ->  data/tetrabonacci_d2_natural_lengths.csv

Figures
-------
figures/tetra_T50_retention.png   — bar chart of retention at T=50
figures/tetra_FSS_lambda1p5.png   — log-log IPR vs N at T=10^4
figures/tetra_D2_natural.png      — log-log IPR vs N at lambda=0 with D₂ fit

Author
------
    Pablo Nogueira Grossi  |  ORCID: 0009-0000-6496-2186
    G6 LLC, Newark NJ  |  pablogrossi@hotmail.com
    GitHub: https://github.com/TOTOGT/AXLE

License: MIT
"""

from __future__ import annotations

import csv
import sys
import time as _time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigh, eigh_tridiagonal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO = Path(__file__).parent
DATA_DIR = REPO / "data"
FIG_DIR = REPO / "figures"
DATA_DIR.mkdir(exist_ok=True)
FIG_DIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Substitution word generators
# ---------------------------------------------------------------------------

TETRA_RULES = {"A": "AB", "B": "AC", "C": "AD", "D": "A"}


def tetrabonacci_word(N: int) -> str:
    """Return the tetrabonacci word truncated to N sites."""
    s = "A"
    while len(s) < N:
        s = "".join(TETRA_RULES[ch] for ch in s)
    return s[:N]


def tetrabonacci_word_natural(n_iterations: int) -> str:
    """Return the n-th iteration of the tetrabonacci substitution from 'A'."""
    s = "A"
    for _ in range(n_iterations):
        s = "".join(TETRA_RULES[ch] for ch in s)
    return s


def fibonacci_word(N: int) -> list[int]:
    """Fibonacci substitution: A->AB, B->A (encoded 0,1)."""
    word = [0]
    rules = {0: [0, 1], 1: [0]}
    while len(word) < N:
        word = [s for c in word for s in rules[c]]
    return word[:N]


def tribonacci_word(N: int) -> list[int]:
    """Tribonacci substitution: 0->01, 1->02, 2->0 (Rauzy)."""
    word = [0]
    rules = {0: [0, 1], 1: [0, 2], 2: [0]}
    while len(word) < N:
        word = [s for c in word for s in rules[c]]
    return word[:N]


# ---------------------------------------------------------------------------
# Hamiltonian builders
# ---------------------------------------------------------------------------

def build_hamiltonian_fib(N: int, t_mod: float = 0.5) -> tuple[np.ndarray, np.ndarray]:
    word = fibonacci_word(N + 1)
    hop_map = {0: 1.0, 1: t_mod}
    hoppings = np.array([hop_map[word[j]] for j in range(N - 1)])
    H = np.zeros((N, N))
    for j in range(N - 1):
        H[j, j + 1] = hoppings[j]
        H[j + 1, j] = hoppings[j]
    return H, hoppings


def build_hamiltonian_trib(N: int, t_mod: float = 0.5) -> tuple[np.ndarray, np.ndarray]:
    word = tribonacci_word(N + 1)
    hop_map = {0: 1.0, 1: t_mod, 2: t_mod ** 2}
    hoppings = np.array([hop_map[word[j]] for j in range(N - 1)])
    H = np.zeros((N, N))
    for j in range(N - 1):
        H[j, j + 1] = hoppings[j]
        H[j + 1, j] = hoppings[j]
    return H, hoppings


def build_hamiltonian_tetra(N: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Tetrabonacci Hamiltonian with hopping convention t_A=1.0, t_B=0.5,
    t_C=0.25, t_D=0.125 (geometric, ratio 0.5 — a convention choice).
    """
    word = tetrabonacci_word(N + 1)
    hop_map = {"A": 1.0, "B": 0.5, "C": 0.25, "D": 0.125}
    hoppings = np.array([hop_map[word[j]] for j in range(N - 1)])
    H = np.zeros((N, N))
    for j in range(N - 1):
        H[j, j + 1] = hoppings[j]
        H[j + 1, j] = hoppings[j]
    return H, hoppings


# ---------------------------------------------------------------------------
# Eigenstate utilities
# ---------------------------------------------------------------------------

def mid_gap_state(H: np.ndarray) -> tuple[np.ndarray, float]:
    """Eigenstate of H closest to E=0."""
    vals, vecs = eigh(H)
    idx = np.argmin(np.abs(vals))
    return vecs[:, idx], float(vals[idx])


def ipr(psi: np.ndarray) -> float:
    """Inverse participation ratio."""
    norm2 = float(np.sum(np.abs(psi) ** 2))
    return float(np.sum(np.abs(psi) ** 4)) / norm2 ** 2


# ---------------------------------------------------------------------------
# DNLS integrators
# ---------------------------------------------------------------------------

def _dnls_rhs(t: float, z: np.ndarray, lam: float, hoppings: np.ndarray) -> np.ndarray:
    N = len(z) >> 1
    x, y = z[:N], z[N:]
    Hx = np.zeros(N)
    Hx[:-1] += hoppings * x[1:]
    Hx[1:] += hoppings * x[:-1]
    Hy = np.zeros(N)
    Hy[:-1] += hoppings * y[1:]
    Hy[1:] += hoppings * y[:-1]
    nl = x * x + y * y
    return np.concatenate([-Hy + lam * nl * y, Hx - lam * nl * x])


def evolve_t50(
    psi0: np.ndarray,
    hoppings: np.ndarray,
    lam: float,
    T: float = 50.0,
    rtol: float = 1e-6,
    atol: float = 1e-8,
) -> tuple[float, float, float]:
    """
    T=50 integration with RK45 (paper-matching tolerances).
    Returns (ipr_t0, ipr_tT, norm_final).
    """
    N = len(psi0)
    psi0 = psi0 / np.sqrt(np.dot(psi0, psi0))
    z0 = np.concatenate([psi0, np.zeros(N)])
    ipr0 = ipr(psi0)

    sol = solve_ivp(
        _dnls_rhs,
        [0.0, T],
        z0,
        method="RK45",
        args=(lam, hoppings),
        rtol=rtol,
        atol=atol,
        max_step=0.1,
        t_eval=[T],
        dense_output=False,
    )
    zf = sol.y[:, -1]
    psi_f = zf[:N] + 1j * zf[N:]
    norm_f = float(np.sqrt(np.sum(np.abs(psi_f) ** 2)))
    ipr_f = ipr(psi_f)
    return ipr0, ipr_f, norm_f


def evolve_long(
    psi0: np.ndarray,
    hoppings: np.ndarray,
    lam: float,
    T: float = 1e4,
    n_checkpoints: int = 1,
    rtol: float = 1e-9,
    atol: float = 1e-11,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Long-time integration with DOP853. Returns (t_arr, ipr_arr, norm_arr).
    If n_checkpoints > 1, uses log-spaced t_eval; otherwise just endpoint.
    """
    N = len(psi0)
    psi0 = psi0 / np.sqrt(np.dot(psi0, psi0))
    z0 = np.concatenate([psi0, np.zeros(N)])
    ipr0 = ipr(psi0)

    if n_checkpoints > 1:
        t_eval = np.geomspace(1.0, T, n_checkpoints)
        t_eval = np.unique(np.concatenate([[0.0], t_eval]))
    else:
        t_eval = np.array([0.0, T])

    sol = solve_ivp(
        _dnls_rhs,
        [0.0, T],
        z0,
        method="DOP853",
        args=(lam, hoppings),
        rtol=rtol,
        atol=atol,
        t_eval=t_eval,
        dense_output=False,
    )

    n_pts = sol.y.shape[1]
    ipr_arr = np.empty(n_pts)
    norm_arr = np.empty(n_pts)
    for k in range(n_pts):
        psi_k = sol.y[:N, k] + 1j * sol.y[N:, k]
        norm_arr[k] = float(np.sqrt(np.sum(np.abs(psi_k) ** 2)))
        ipr_arr[k] = ipr(psi_k)

    # Overwrite t=0 slot with the analytically computed initial-state values.
    # solve_ivp evaluates the RHS at t=0 during setup, which can introduce
    # tiny numerical artefacts before any actual time stepping has occurred.
    ipr_arr[0] = ipr0
    norm_arr[0] = 1.0

    return sol.t, ipr_arr, norm_arr


# ---------------------------------------------------------------------------
# Task B — OEIS A000078 verification
# ---------------------------------------------------------------------------

def task_b_verify() -> None:
    print("=" * 70)
    print("[1] OEIS A000078 VERIFICATION")
    print("=" * 70)
    expected = [1, 2, 4, 8, 15, 29, 56, 108, 208, 401, 773, 1490, 2872, 5536]
    all_ok = True
    for n in range(14):
        s = tetrabonacci_word_natural(n)
        L = len(s)
        e = expected[n]
        status = "OK" if L == e else "FAIL"
        if L != e:
            all_ok = False
        print(f"  tetrabonacci_word_natural({n:>2})  length = {L:<6} (expected {e:<6}) {status}")
    if not all_ok:
        print("\nHALT: length mismatch against OEIS A000078 — halting for human review.")
        sys.exit(1)
    print("\n  All 14 lengths verified against OEIS A000078. Proceeding.\n")


# ---------------------------------------------------------------------------
# Task C — T=50, N=500
# ---------------------------------------------------------------------------

LAMBDAS_C = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0]
N_C = 500
T_C = 50.0
NORM_TOL = 1e-5


def task_c_t50() -> dict:
    """
    Run T=50 for fibonacci, tribonacci, tetrabonacci at N=500.
    Returns retention dict keyed by (chain, lambda).
    """
    print("=" * 70)
    print("[2] TASK C — T=50 retention at N=500")
    print("=" * 70)

    # Build chains
    H_fib, hop_fib = build_hamiltonian_fib(N_C)
    H_trib, hop_trib = build_hamiltonian_trib(N_C)
    H_tetra, hop_tetra = build_hamiltonian_tetra(N_C)

    psi0_fib, E0_fib = mid_gap_state(H_fib)
    psi0_trib, E0_trib = mid_gap_state(H_trib)
    psi0_tetra, E0_tetra = mid_gap_state(H_tetra)

    print(f"  fibonacci  mid-gap E0={E0_fib:.6f}  IPR(0)={ipr(psi0_fib):.6f}")
    print(f"  tribonacci mid-gap E0={E0_trib:.6f}  IPR(0)={ipr(psi0_trib):.6f}")
    print(f"  tetra      mid-gap E0={E0_tetra:.6f}  IPR(0)={ipr(psi0_tetra):.6f}")
    print()

    chains = [
        ("fibonacci", psi0_fib, hop_fib),
        ("tribonacci", psi0_trib, hop_trib),
        ("tetrabonacci", psi0_tetra, hop_tetra),
    ]

    rows: list[dict] = []
    retention: dict[tuple[str, float], float] = {}
    max_norm_leak = 0.0
    norm_flags = []

    n_total = len(chains) * len(LAMBDAS_C)
    run_idx = 0
    for chain_name, psi0, hoppings in chains:
        for lam in LAMBDAS_C:
            run_idx += 1
            t0 = _time.perf_counter()
            ipr0, ipr_T, norm_f = evolve_t50(psi0, hoppings, lam, T=T_C)
            elapsed = _time.perf_counter() - t0
            ret = ipr_T / ipr0
            retention[(chain_name, lam)] = ret
            leak = abs(norm_f - 1.0)
            max_norm_leak = max(max_norm_leak, leak)
            flag = ""
            if leak > NORM_TOL:
                flag = "  *** NORM LEAK ***"
                norm_flags.append((chain_name, lam, leak))
            print(
                f"  [{run_idx:>2}/{n_total}] {chain_name:<14} lambda={lam:<5.1f} "
                f"IPR(0)={ipr0:.6f}  IPR(T)={ipr_T:.6f}  retention={ret:.4f}"
                f"  norm_leak={leak:.2e}  ({elapsed:.1f}s){flag}"
            )
            rows.append({
                "time": T_C,
                "lambda": lam,
                "chain": chain_name,
                "IPR": ipr_T,
                "norm": norm_f,
            })

    # Write CSV
    csv_path = DATA_DIR / "tetrabonacci_T50_N500.csv"
    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["time", "lambda", "chain", "IPR", "norm"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"\n  Wrote {len(rows)} rows -> {csv_path}")
    print(f"  Max norm leak in T=50 runs: {max_norm_leak:.2e}")
    if norm_flags:
        for ch, lam, lk in norm_flags:
            print(f"  FLAG: {ch} lambda={lam} leak={lk:.2e}")
    else:
        print("  Flags raised: none")

    # Print retention table
    print()
    print("  Retention table (IPR(T=50) / IPR(t=0)):")
    print(f"  {'lambda':>6}  {'fib':>12}  {'trib':>12}  {'tetra':>12}")
    print("  " + "-" * 48)
    for lam in LAMBDAS_C:
        r_fib = retention.get(("fibonacci", lam), float("nan"))
        r_trib = retention.get(("tribonacci", lam), float("nan"))
        r_tetra = retention.get(("tetrabonacci", lam), float("nan"))
        print(f"  {lam:>6.1f}  {r_fib:>12.6f}  {r_trib:>12.6f}  {r_tetra:>12.6f}")

    return {"retention": retention, "max_norm_leak": max_norm_leak, "norm_flags": norm_flags}


# ---------------------------------------------------------------------------
# Task D — T=10^4 FSS at lambda=1.5
# ---------------------------------------------------------------------------

N_D_SIZES = [500, 1000, 2000]
LAM_D = 1.5
T_D = 1e4


def task_d_fss() -> dict:
    print()
    print("=" * 70)
    print("[3] TASK D — T=10^4 FSS at lambda=1.5")
    print("=" * 70)

    rows: list[dict] = []
    ipr_results: dict[tuple[str, int], float] = {}
    max_norm_leak = 0.0
    norm_flags = []

    chain_builders = [
        ("fibonacci", build_hamiltonian_fib),
        ("tribonacci", build_hamiltonian_trib),
        ("tetrabonacci", lambda N: build_hamiltonian_tetra(N)),
    ]

    n_total = len(chain_builders) * len(N_D_SIZES)
    run_idx = 0
    for chain_name, builder in chain_builders:
        for N in N_D_SIZES:
            run_idx += 1
            H, hoppings = builder(N)
            psi0, E0 = mid_gap_state(H)
            ipr0 = ipr(psi0)
            t0 = _time.perf_counter()
            t_arr, ipr_arr, norm_arr = evolve_long(psi0, hoppings, LAM_D, T=T_D)
            elapsed = _time.perf_counter() - t0
            ipr_T = float(ipr_arr[-1])
            norm_T = float(norm_arr[-1])
            leak = float(np.max(np.abs(norm_arr - 1.0)))
            max_norm_leak = max(max_norm_leak, leak)
            ipr_results[(chain_name, N)] = ipr_T
            flag = ""
            if leak > NORM_TOL:
                flag = "  *** NORM LEAK ***"
                norm_flags.append((chain_name, N, leak))
            print(
                f"  [{run_idx:>2}/{n_total}] {chain_name:<14} N={N:<5} "
                f"IPR(0)={ipr0:.6f}  IPR(T)={ipr_T:.6f}  leak={leak:.2e}"
                f"  ({elapsed:.1f}s){flag}"
            )
            for t_k, ipr_k, norm_k in zip(t_arr, ipr_arr, norm_arr):
                rows.append({
                    "time": t_k,
                    "lambda": LAM_D,
                    "chain": chain_name,
                    "N": N,
                    "IPR": ipr_k,
                    "norm": norm_k,
                })

    csv_path = DATA_DIR / "tetrabonacci_lambda1p5_T1e4.csv"
    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["time", "lambda", "chain", "N", "IPR", "norm"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"\n  Wrote {len(rows)} rows -> {csv_path}")
    print(f"  Max norm leak in T=10^4 runs: {max_norm_leak:.2e}")
    if norm_flags:
        for ch, N, lk in norm_flags:
            print(f"  FLAG: {ch} N={N} leak={lk:.2e}")
    else:
        print("  Flags raised: none")

    # Print FSS table
    print()
    print("  T=10^4 FSS table at lambda=1.5:")
    print(f"  {'':>16} {'N=500':>12}  {'N=1000':>12}  {'N=2000':>12}")
    print("  " + "-" * 56)
    for ch in ["fibonacci", "tribonacci", "tetrabonacci"]:
        vals = [ipr_results.get((ch, N), float("nan")) for N in N_D_SIZES]
        print(f"  {ch + ' IPR':<16} {vals[0]:>12.6f}  {vals[1]:>12.6f}  {vals[2]:>12.6f}")
    # Ratio rows
    for (ca, cb) in [("fibonacci", "tribonacci"), ("tribonacci", "tetrabonacci"), ("fibonacci", "tetrabonacci")]:
        label = f"{ca[:4]}/{cb[:4]}"
        vals = [ipr_results.get((ca, N), float("nan")) / ipr_results.get((cb, N), float("nan"))
                for N in N_D_SIZES]
        print(f"  {label:<16} {vals[0]:>12.6f}  {vals[1]:>12.6f}  {vals[2]:>12.6f}")

    return {
        "ipr_results": ipr_results,
        "rows": rows,
        "max_norm_leak": max_norm_leak,
        "norm_flags": norm_flags,
    }


# ---------------------------------------------------------------------------
# Task E — D₂ from mid-gap IPR at natural lengths, lambda=0
# ---------------------------------------------------------------------------

# Natural lengths correspond to iterations 8..14 (OEIS A000078)
NATURAL_ITERS_E = [8, 9, 10, 11, 12, 13, 14]
NATURAL_LENGTHS_E = [208, 401, 773, 1490, 2872, 5536, 10671]
SPREAD_FRAC = 0.03  # spatial spread filter: sigma >= SPREAD_FRAC * N


def _spatial_spread(psi: np.ndarray) -> float:
    """Return the spatial standard deviation of |psi|^2."""
    prob = np.abs(psi) ** 2
    prob /= prob.sum()
    j = np.arange(len(psi), dtype=float)
    mu = float(np.dot(prob, j))
    return float(np.sqrt(np.dot(prob, (j - mu) ** 2)))


def task_e_d2() -> dict:
    print()
    print("=" * 70)
    print("[4] TASK E — D₂ from mid-gap eigenstate IPR at natural lengths")
    print("=" * 70)

    rows_tetra: list[dict] = []
    n_filtered_total = 0

    # Also compute fib / trib for comparison
    rows_fib: list[dict] = []
    rows_trib: list[dict] = []

    # Fibonacci natural lengths: Fibonacci numbers >= 208 within range
    fib_lengths = [233, 377, 610, 987, 1597, 2584]
    trib_lengths = [274, 504, 927, 1705, 3136]

    print("\n  Tetrabonacci natural lengths:")
    print(f"  {'chain':<14} {'N':>6}  {'IPR':>12}  {'log10(N)':>10}  {'log10(IPR)':>12}  {'filtered':>8}")
    print("  " + "-" * 68)

    for it, N in zip(NATURAL_ITERS_E, NATURAL_LENGTHS_E):
        word_str = tetrabonacci_word_natural(it)
        hop_map = {"A": 1.0, "B": 0.5, "C": 0.25, "D": 0.125}
        # use eigh_tridiagonal for efficiency
        diag = np.zeros(N)
        off = np.array([hop_map[word_str[j]] for j in range(N - 1)])
        vals, vecs = eigh_tridiagonal(diag, off)
        # find mid-gap state(s), apply spread filter
        order = np.argsort(np.abs(vals))
        n_filtered = 0
        selected_ipr = None
        for idx in order:
            psi = vecs[:, idx]
            sigma = _spatial_spread(psi)
            if sigma < SPREAD_FRAC * N:
                # Exclude compact zero-modes: states localised on a few
                # sites do not represent the extended multifractal behaviour
                # and would corrupt the D₂ fit.
                n_filtered += 1
                continue
            selected_ipr = ipr(psi)
            break
        n_filtered_total += n_filtered
        if selected_ipr is None:
            print(f"  {'tetrabonacci':<14} {N:>6}  {'NO STATE PASSED FILTER':>12}")
            continue
        log_N = float(np.log10(N))
        log_ipr = float(np.log10(selected_ipr))
        print(f"  {'tetrabonacci':<14} {N:>6}  {selected_ipr:>12.8f}  {log_N:>10.4f}  {log_ipr:>12.6f}  {n_filtered:>8d}")
        rows_tetra.append({
            "chain": "tetrabonacci",
            "N": N,
            "IPR": selected_ipr,
            "log10_N": log_N,
            "log10_IPR": log_ipr,
        })

    # Fit D₂ for tetrabonacci
    if len(rows_tetra) >= 2:
        x = np.array([r["log10_N"] for r in rows_tetra])
        y = np.array([r["log10_IPR"] for r in rows_tetra])
        p, cov = np.polyfit(x, y, 1, cov=True)
        d2_tetra = float(-p[0])
        d2_tetra_err = float(np.sqrt(cov[0, 0]))
        y_pred = np.polyval(p, x)
        ss_res = float(np.sum((y - y_pred) ** 2))
        ss_tot = float(np.sum((y - y.mean()) ** 2))
        r2_tetra = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
        n_pts_tetra = len(rows_tetra)
    else:
        d2_tetra = d2_tetra_err = r2_tetra = float("nan")
        n_pts_tetra = len(rows_tetra)

    # Fibonacci reference
    print("\n  Fibonacci reference:")
    print(f"  {'chain':<14} {'N':>6}  {'IPR':>12}  {'log10(N)':>10}  {'log10(IPR)':>12}")
    print("  " + "-" * 60)
    for N in fib_lengths:
        word = fibonacci_word(N + 1)
        hop_map_fib = {0: 1.0, 1: 0.5}
        diag = np.zeros(N)
        off = np.array([hop_map_fib[word[j]] for j in range(N - 1)])
        vals, vecs = eigh_tridiagonal(diag, off)
        order = np.argsort(np.abs(vals))
        for idx in order:
            psi = vecs[:, idx]
            if _spatial_spread(psi) >= SPREAD_FRAC * N:
                v = ipr(psi)
                log_N = float(np.log10(N))
                log_ipr = float(np.log10(v))
                print(f"  {'fibonacci':<14} {N:>6}  {v:>12.8f}  {log_N:>10.4f}  {log_ipr:>12.6f}")
                rows_fib.append({"chain": "fibonacci", "N": N, "IPR": v,
                                 "log10_N": log_N, "log10_IPR": log_ipr})
                break

    # Tribonacci reference
    print("\n  Tribonacci reference:")
    print(f"  {'chain':<14} {'N':>6}  {'IPR':>12}  {'log10(N)':>10}  {'log10(IPR)':>12}")
    print("  " + "-" * 60)
    for N in trib_lengths:
        word = tribonacci_word(N + 1)
        hop_map_trib = {0: 1.0, 1: 0.5, 2: 0.25}
        diag = np.zeros(N)
        off = np.array([hop_map_trib[word[j]] for j in range(N - 1)])
        vals, vecs = eigh_tridiagonal(diag, off)
        order = np.argsort(np.abs(vals))
        for idx in order:
            psi = vecs[:, idx]
            if _spatial_spread(psi) >= SPREAD_FRAC * N:
                v = ipr(psi)
                log_N = float(np.log10(N))
                log_ipr = float(np.log10(v))
                print(f"  {'tribonacci':<14} {N:>6}  {v:>12.8f}  {log_N:>10.4f}  {log_ipr:>12.6f}")
                rows_trib.append({"chain": "tribonacci", "N": N, "IPR": v,
                                  "log10_N": log_N, "log10_IPR": log_ipr})
                break

    # Fit D₂ for fib and trib
    def fit_d2(rows):
        if len(rows) < 2:
            return float("nan"), float("nan"), float("nan"), len(rows)
        x = np.array([r["log10_N"] for r in rows])
        y = np.array([r["log10_IPR"] for r in rows])
        p, cov = np.polyfit(x, y, 1, cov=True)
        d2 = float(-p[0])
        d2_err = float(np.sqrt(cov[0, 0]))
        y_pred = np.polyval(p, x)
        ss_res = float(np.sum((y - y_pred) ** 2))
        ss_tot = float(np.sum((y - y.mean()) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
        return d2, d2_err, r2, len(rows)

    d2_fib, d2_fib_err, r2_fib, n_fib = fit_d2(rows_fib)
    d2_trib, d2_trib_err, r2_trib, n_trib = fit_d2(rows_trib)

    # Write CSV (tetrabonacci only per brief)
    all_rows = rows_tetra + rows_fib + rows_trib
    csv_path = DATA_DIR / "tetrabonacci_d2_natural_lengths.csv"
    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["chain", "N", "IPR", "log10_N", "log10_IPR"])
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"\n  Wrote {len(all_rows)} rows -> {csv_path}")
    print(f"  Compact zero-modes filtered: {n_filtered_total} total (spread < {SPREAD_FRAC}*N criterion)")

    print()
    print("  D₂ fits:")
    print(f"  D₂_fib   = {d2_fib:.4f}  std_err = {d2_fib_err:.4f}  R² = {r2_fib:.4f}  n_pts = {n_fib}")
    print(f"  D₂_trib  = {d2_trib:.4f}  std_err = {d2_trib_err:.4f}  R² = {r2_trib:.4f}  n_pts = {n_trib}")
    print(f"  D₂_tetra = {d2_tetra:.4f}  std_err = {d2_tetra_err:.4f}  R² = {r2_tetra:.4f}  n_pts = {n_pts_tetra}")
    print()
    print(f"  Comparison: D₂_fib = {d2_fib:.3f}, D₂_trib = {d2_trib:.3f}, D₂_tetra = {d2_tetra:.3f}")

    return {
        "rows_tetra": rows_tetra,
        "rows_fib": rows_fib,
        "rows_trib": rows_trib,
        "d2_fib": d2_fib, "d2_fib_err": d2_fib_err, "r2_fib": r2_fib, "n_fib": n_fib,
        "d2_trib": d2_trib, "d2_trib_err": d2_trib_err, "r2_trib": r2_trib, "n_trib": n_trib,
        "d2_tetra": d2_tetra, "d2_tetra_err": d2_tetra_err, "r2_tetra": r2_tetra, "n_pts_tetra": n_pts_tetra,
        "n_filtered_total": n_filtered_total,
    }


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def make_figures(c_ret: dict, d_res: dict, e_res: dict) -> None:
    print()
    print("=" * 70)
    print("[5] GENERATING FIGURES")
    print("=" * 70)

    retention = c_ret["retention"]
    ipr_fss = d_res["ipr_results"]

    # --- Figure 1: T=50 retention bar chart ---
    fig_path = FIG_DIR / "tetra_T50_retention.png"
    chains_order = ["fibonacci", "tribonacci", "tetrabonacci"]
    chain_labels = ["Fib", "Trib", "Tetra"]
    colors = ["#4878CF", "#6ACC65", "#D65F5F"]
    x = np.arange(len(LAMBDAS_C))
    width = 0.25

    fig, ax = plt.subplots(figsize=(12, 5), dpi=140)
    for i, (ch, label, col) in enumerate(zip(chains_order, chain_labels, colors)):
        vals = [retention.get((ch, lam), float("nan")) for lam in LAMBDAS_C]
        ax.bar(x + i * width, vals, width, label=label, color=col, alpha=0.85)
    ax.set_xlabel("λ (nonlinearity)")
    ax.set_ylabel("Retention = IPR(T=50) / IPR(0)")
    ax.set_title("T=50 IPR Retention: Fibonacci vs Tribonacci vs Tetrabonacci (N=500)")
    ax.set_xticks(x + width)
    ax.set_xticklabels([str(l) for l in LAMBDAS_C])
    ax.legend()
    ax.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(fig_path, dpi=140)
    plt.close()
    print(f"  -> {fig_path}")

    # --- Figure 2: T=10^4 FSS log-log ---
    fig_path = FIG_DIR / "tetra_FSS_lambda1p5.png"
    fig, ax = plt.subplots(figsize=(7, 5), dpi=140)
    markers = {"fibonacci": "o", "tribonacci": "s", "tetrabonacci": "^"}
    colors2 = {"fibonacci": "#4878CF", "tribonacci": "#6ACC65", "tetrabonacci": "#D65F5F"}
    for ch in chains_order:
        Ns = N_D_SIZES
        iprs = [ipr_fss.get((ch, N), float("nan")) for N in Ns]
        ax.loglog(Ns, iprs, marker=markers[ch], color=colors2[ch],
                  lw=2, ms=8, label=ch.capitalize())
    ax.set_xlabel("N (chain length)")
    ax.set_ylabel("IPR(T=10⁴)")
    ax.set_title("FSS: IPR at T=10⁴, λ=1.5")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    plt.savefig(fig_path, dpi=140)
    plt.close()
    print(f"  -> {fig_path}")

    # --- Figure 3: D₂ log-log with fit lines ---
    fig_path = FIG_DIR / "tetra_D2_natural.png"
    fig, ax = plt.subplots(figsize=(7, 5), dpi=140)

    def plot_d2_chain(rows, d2, label, color, marker):
        if not rows:
            return
        Ns = np.array([r["N"] for r in rows])
        iprs = np.array([r["IPR"] for r in rows])
        ax.loglog(Ns, iprs, marker=marker, color=color, lw=0, ms=8,
                  label=f"{label} data")
        if np.isfinite(d2):
            # fit line
            log_Ns = np.log10(Ns)
            log_iprs = np.array([r["log10_IPR"] for r in rows])
            intercept = np.mean(log_iprs) + d2 * np.mean(log_Ns)
            N_fit = np.array([Ns[0] * 0.8, Ns[-1] * 1.2])
            ipr_fit = 10 ** (intercept - d2 * np.log10(N_fit))
            ax.loglog(N_fit, ipr_fit, "--", color=color, lw=2,
                      label=f"{label} D₂={d2:.3f}")

    plot_d2_chain(e_res["rows_fib"], e_res["d2_fib"], "Fib", "#4878CF", "o")
    plot_d2_chain(e_res["rows_trib"], e_res["d2_trib"], "Trib", "#6ACC65", "s")
    plot_d2_chain(e_res["rows_tetra"], e_res["d2_tetra"], "Tetra", "#D65F5F", "^")

    ax.set_xlabel("N (chain length)")
    ax.set_ylabel("IPR (λ=0, mid-gap eigenstate)")
    ax.set_title("D₂ multifractal dimension: IPR ~ N^(-D₂)")
    ax.legend(fontsize=9)
    ax.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    plt.savefig(fig_path, dpi=140)
    plt.close()
    print(f"  -> {fig_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    print("\ndnls_tetrabonacci.py — Tetrabonacci extension")
    print("=" * 70)

    # Task B
    task_b_verify()

    # Task C
    c_ret = task_c_t50()

    # Task D
    d_res = task_d_fss()

    # Task E
    e_res = task_e_d2()

    # Figures
    make_figures(c_ret, d_res, e_res)

    # Summary
    print()
    print("=" * 70)
    print("[6] NORM CONSERVATION SUMMARY")
    print("=" * 70)
    print(f"  Max leak in T=50 runs:   {c_ret['max_norm_leak']:.2e}")
    print(f"  Max leak in T=10^4 runs: {d_res['max_norm_leak']:.2e}")
    t50_flags = len(c_ret['norm_flags'])
    t1e4_flags = len(d_res['norm_flags'])
    print(f"  Flags raised (T=50):     {t50_flags}")
    print(f"  Flags raised (T=10^4):   {t1e4_flags}")

    print()
    print("=" * 70)
    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
