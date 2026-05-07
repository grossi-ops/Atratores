#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import sys
import time
import warnings
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eig
from scipy.optimize import curve_fit
from scipy.optimize import OptimizeWarning

from DNLS.dnls_nbonacci import nbonacci_word, pentabonacci_word_natural
from DNLS.diffusion_solver import (
    build_fission_matrix,
    build_loss_matrix,
    keff,
    lambda_c,
    material_params,
    run_smoke_test_or_fail,
)

DATA_DIR = REPO_ROOT / "data"
FIG_DIR = REPO_ROOT / "figures"
OUT_CSV = DATA_DIR / "lambdac_sweep.csv"
OUT_REPORT = DATA_DIR / "criticality_report.txt"

COL_FIB = "#2166ac"
COL_TRIB = "#d6604d"
COL_GREEN = "#4dac26"
COL_GOLD = "#c9a84c"

GRID = {
    2: [7, 8, 9, 10],
    3: [6, 7, 8, 9],
    4: [6, 7, 8],
    5: [6, 7, 8],
}

COL_BY_N = {2: COL_FIB, 3: COL_TRIB, 4: COL_GREEN, 5: COL_GOLD}


def ensure_dirs() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)


def verify_pentabonacci_oeis() -> list[str]:
    expected = [1, 2, 4, 8, 16, 31, 61, 120, 236]
    lines: list[str] = []
    for i, exp in enumerate(expected):
        got = len(pentabonacci_word_natural(i))
        status = "OK" if got == exp else "FAIL"
        line = (
            f"pentabonacci_word_natural({i})  length = {got:<5d} "
            f"(expected {exp})  {status}"
        )
        print(line)
        lines.append(line)
        if got != exp:
            raise RuntimeError("Pentabonacci OEIS A001591 mismatch; halting.")
    return lines


def spectral_gap(n: int) -> tuple[float, float, float]:
    coeffs = [1] + [-1] * n
    roots = np.roots(coeffs)
    mags = np.sort(np.abs(roots))[::-1]
    rho = float(mags[0])
    rho2 = float(mags[1])
    return rho, rho2, rho - rho2


def dominant_mode(L: np.ndarray, F: np.ndarray) -> np.ndarray:
    vals, vecs = eig(F, L)
    mask = np.isfinite(vals) & (np.abs(vals.imag) < 1e-8)
    vals_r = vals[mask].real
    vecs_r = vecs[:, mask]
    idx = int(np.argmax(vals_r))
    mode = np.real(vecs_r[:, idx])
    if np.max(np.abs(mode)) == 0:
        return mode
    if np.sum(mode) < 0:
        mode = -mode
    return mode / np.max(np.abs(mode))


def exp_model(g: np.ndarray, l_inf: float, C: float, tau: float) -> np.ndarray:
    return l_inf + C * np.exp(-g / tau)


def run() -> dict:
    ensure_dirs()
    smoke = run_smoke_test_or_fail()
    oeis_lines = verify_pentabonacci_oeis()

    sweep_rows: list[dict] = []
    lam_by_n: dict[int, list[tuple[int, int, float]]] = {n: [] for n in GRID}
    mode_by_n: dict[int, np.ndarray] = {}

    for n, gs in GRID.items():
        for g in gs:
            word = nbonacci_word(n, g)
            N = len(word)
            t0 = time.perf_counter()
            lc = lambda_c(word, bracket=(0.3, 6.0), tol=1e-9)
            D, Sigma_r, nuSigmaf = material_params(word, lc)
            L = build_loss_matrix(D, Sigma_r, h=1.0)
            F = build_fission_matrix(nuSigmaf)
            k = keff(L, F)
            wall = time.perf_counter() - t0
            sweep_rows.append(
                {
                    "n": n,
                    "g": g,
                    "N": N,
                    "lambdac": lc,
                    "keff_at_lambdac": k,
                    "wallclock_s": wall,
                }
            )
            lam_by_n[n].append((g, N, lc))

    for n in GRID:
        g_max = max(GRID[n])
        word = nbonacci_word(n, g_max)
        lc = [x[2] for x in lam_by_n[n] if x[0] == g_max][0]
        D, Sigma_r, nuSigmaf = material_params(word, lc)
        L = build_loss_matrix(D, Sigma_r, h=1.0)
        F = build_fission_matrix(nuSigmaf)
        mode_by_n[n] = dominant_mode(L, F)

    with OUT_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["n", "g", "N", "lambdac", "keff_at_lambdac", "wallclock_s"],
        )
        writer.writeheader()
        for row in sorted(sweep_rows, key=lambda r: (r["n"], r["g"])):
            writer.writerow(row)

    gaps = {n: spectral_gap(n) for n in GRID}
    lam_lim = {n: sorted(vals, key=lambda x: x[0])[-1][2] for n, vals in lam_by_n.items()}
    deltas = np.array([gaps[n][2] for n in [2, 3, 4, 5]], dtype=float)
    lambdas = np.array([lam_lim[n] for n in [2, 3, 4, 5]], dtype=float)

    fit_params, fit_cov = np.polyfit(deltas, lambdas, 1, cov=True)
    alpha, beta = fit_params
    alpha_err = float(np.sqrt(fit_cov[0, 0]))
    beta_err = float(np.sqrt(fit_cov[1, 1]))
    pred = alpha * deltas + beta
    residuals = lambdas - pred
    r = float(np.corrcoef(deltas, lambdas)[0, 1])
    r2 = r * r

    tau_fit: dict[int, float] = {}
    tau_pred: dict[int, float] = {}
    for n, vals in lam_by_n.items():
        arr = np.array(sorted(vals, key=lambda x: x[0]), dtype=float)
        g = arr[:, 0]
        y = arr[:, 2]
        if float(np.max(np.abs(y - y[-1]))) < 1e-12:
            tau_fit[n] = float("inf")
        else:
            p0 = [y[-1], y[0] - y[-1], 1.0]
            bounds = ([0.0, -10.0, 1e-6], [10.0, 10.0, 1e3])
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", OptimizeWarning)
                    popt, _ = curve_fit(exp_model, g, y, p0=p0, bounds=bounds, maxfev=10000)
                tau_fit[n] = float(popt[2])
            except Exception:
                tau_fit[n] = float("nan")
        rho, rho2, _ = gaps[n]
        tau_pred[n] = 1.0 / math.log(rho / rho2)

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["DejaVu Serif", "Times New Roman", "Georgia"],
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 150,
            "savefig.dpi": 300,
        }
    )

    fig1, ax1 = plt.subplots(figsize=(6.0, 4.5))
    xfit = np.linspace(deltas.min() * 0.98, deltas.max() * 1.02, 200)
    yfit = alpha * xfit + beta
    ax1.scatter(deltas, lambdas, s=60, c=[COL_BY_N[n] for n in [2, 3, 4, 5]])
    ax1.plot(xfit, yfit, color="black", lw=1.5, label=f"λ_c = {alpha:.3f}Δ + {beta:.3f}")
    for i, n in enumerate([2, 3, 4, 5]):
        ax1.annotate(f"n={n}", (deltas[i], lambdas[i]), textcoords="offset points", xytext=(5, 5))
    ax1.set_xlabel("Spectral gap Δ_n")
    ax1.set_ylabel("Converged criticality λ_c(n)")
    ax1.legend(frameon=True)
    ax1.grid(alpha=0.3)
    fig1.tight_layout()
    fig1.savefig(FIG_DIR / "fig_criticality_correlation.pdf")
    plt.close(fig1)

    fig2, ax2 = plt.subplots(figsize=(7.0, 4.8))
    for n in [2, 3, 4, 5]:
        arr = np.array(sorted(lam_by_n[n], key=lambda x: x[0]), dtype=float)
        ax2.plot(arr[:, 0], arr[:, 2], "o-", color=COL_BY_N[n], label=f"n={n}")
    ax2.axhline(1.064, color=COL_FIB, ls="--", lw=1.0, alpha=0.6, label="target n=2: 1.064")
    ax2.axhline(37.0 / 32.0, color=COL_TRIB, ls="--", lw=1.0, alpha=0.6, label="target n=3: 37/32")
    ax2.axhline(7.0 / 6.0, color=COL_GREEN, ls="--", lw=1.0, alpha=0.6, label="target n≥4: 7/6")
    if 5 in lam_lim:
        ax2.axhline(
            lam_lim[5],
            color=COL_GOLD,
            ls=":",
            lw=1.0,
            alpha=0.8,
            label=f"n=5 limit: {lam_lim[5]:.6f}",
        )
    ax2.set_xlabel("Generation g")
    ax2.set_ylabel("λ_c(n,g)")
    ax2.grid(alpha=0.3)
    ax2.legend(ncol=2, fontsize=8, frameon=True)
    fig2.tight_layout()
    fig2.savefig(FIG_DIR / "fig_criticality_saturation.pdf")
    plt.close(fig2)

    fig3, axs = plt.subplots(2, 2, figsize=(8.0, 5.5), sharex=False, sharey=True)
    for ax, n in zip(axs.ravel(), [2, 3, 4, 5]):
        phi = mode_by_n[n]
        x = np.arange(phi.size)
        ax.plot(x, phi, color=COL_BY_N[n], lw=1.2)
        ax.set_title(f"n={n}, g_max={max(GRID[n])}, N={phi.size}")
        ax.set_xlabel("cell index")
        ax.set_ylabel("φ₀ / max|φ₀|")
        ax.grid(alpha=0.25)
    fig3.tight_layout()
    fig3.savefig(FIG_DIR / "fig_criticality_flux_modes.pdf")
    plt.close(fig3)

    section_2_rows = []
    draft_values = {2: "1.064", 3: "37/32 = 1.15625", 4: "7/6 ≈ 1.16667", 5: "7/6 ≈ 1.16667"}
    gmaxs = {2: 10, 3: 9, 4: 8, 5: 8}
    for n in [2, 3, 4, 5]:
        rho, rho2, delta = gaps[n]
        gmax = gmaxs[n]
        N = int([x[1] for x in lam_by_n[n] if x[0] == gmax][0])
        lc = [x[2] for x in lam_by_n[n] if x[0] == gmax][0]
        section_2_rows.append((n, rho, rho2, delta, gmax, N, lc, draft_values[n]))

    with OUT_REPORT.open("w") as f:
        f.write("### [1] Pentabonacci OEIS verification\n")
        f.write("```\n")
        for line in oeis_lines:
            f.write(line + "\n")
        f.write("```\n\n")

        f.write("### [2] λ_c(n) headline table\n")
        f.write("```\n")
        f.write("n   ρ_n         |ρ_n^(2)|    Δ_n         g_max  N      λ_c(n,g_max)   draft_value\n")
        for row in section_2_rows:
            n, rho, rho2, delta, gmax, N, lc, dval = row
            f.write(f"{n:<1d}   {rho:.5f}     {rho2:.5f}      {delta:.5f}     {gmax:<2d}     {N:<4d}   {lc:.12f}   {dval}\n")
        f.write("```\n\n")

        f.write("### [3] Linear fit\n")
        f.write("```\n")
        f.write("Fit  λ_c = α·Δ + β  over (n=2,3,4,5):\n")
        f.write(f"  α       = {alpha:.12f} ± {alpha_err:.12f}\n")
        f.write(f"  β       = {beta:.12f} ± {beta_err:.12f}\n")
        f.write(f"  r       = {r:.12f}\n")
        f.write(f"  r²      = {r2:.12f}\n")
        f.write(
            "  residuals (per n): "
            f"n=2: {residuals[0]:.12e}, n=3: {residuals[1]:.12e}, "
            f"n=4: {residuals[2]:.12e}, n=5: {residuals[3]:.12e}\n"
        )
        f.write("```\n\n")

        f.write("### [4] Generation convergence\n")
        f.write("```\n")
        for n in [2, 3, 4, 5]:
            arr = sorted(lam_by_n[n], key=lambda x: x[0])
            joined = ", ".join([f"g={int(g)}: {lc:.12f}" for g, _, lc in arr])
            f.write(f"n={n}: {joined}\n")
        f.write("\n")
        f.write("Fitted convergence time τ_n (from λ_c(n,g) = λ_c(n) + C·exp(−g/τ_n)):\n")
        for n in [2, 3, 4, 5]:
            f.write(
                f"n={n}: τ = {tau_fit[n]:.12f}, "
                f"predicted 1/log(ρ_{n}/|ρ_{n}^(2)|) = {tau_pred[n]:.12f}\n"
            )
        f.write("```\n\n")

        f.write("### [5] Smoke test\n")
        f.write(
            f"Uniform-fissile slab (all A_0, N=50): computed λ_c = {smoke.computed_lambda_c:.12f}, "
            f"analytic prediction λ_c = D·(π/L)² + Σ_r = {smoke.analytic_lambda_c:.12f}, "
            f"relative error {smoke.relative_error:.12%}.\n\n"
        )

        f.write("### [6] Three figures committed\n")
        f.write("Paths under `figures/`.\n\n")

        f.write("### [7] Plain-language summary\n")
        f.write(
            "The spectral-gap correlation is positive and approximately linear across n=2..5, "
            "with fit metrics reported above. "
            "The computed λ_c values should be compared directly against 1.064, 37/32, and 7/6 "
            "without forcing rational rounding. "
            "Generation-wise convergence is quantified by fitted τ_n values and compared to "
            "1/log(ρ_n/|ρ_n^(2)|) predictions for each n.\n"
        )

    return {
        "smoke": smoke,
        "oeis_lines": oeis_lines,
        "gaps": gaps,
        "lam_by_n": lam_by_n,
        "lam_lim": lam_lim,
        "fit": {
            "alpha": alpha,
            "beta": beta,
            "alpha_err": alpha_err,
            "beta_err": beta_err,
            "r": r,
            "r2": r2,
            "residuals": residuals,
        },
        "tau_fit": tau_fit,
        "tau_pred": tau_pred,
    }


if __name__ == "__main__":
    out = run()
    print(f"Wrote: {OUT_CSV}")
    print(f"Wrote: {OUT_REPORT}")
    print(f"Wrote figures under: {FIG_DIR}")
    print("Converged λ_c values:", {k: v for k, v in out["lam_lim"].items()})
