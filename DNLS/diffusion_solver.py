#!/usr/bin/env python3
"""
One-group diffusion criticality solver for n-bonacci media.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.linalg import eig
from scipy.optimize import brentq
from scipy.sparse import csc_matrix, diags
from scipy.sparse.linalg import eigs


def mesh_subdivision_factor(h: float) -> int:
    """Return integer m such that h = 1/m, validating exact unit-cell subdivision."""
    if h <= 0:
        raise ValueError("h must be positive")
    inv = 1.0 / h
    m = int(round(inv))
    if m < 1 or not np.isclose(inv, float(m), rtol=0.0, atol=1e-12):
        raise ValueError("mesh h must divide unit cells exactly (use h = 1/k for integer k)")
    return m


def refined_word_for_mesh(word: str, h: float) -> str:
    """Repeat each symbol m=1/h times so each original unit cell is subdivided uniformly."""
    m = mesh_subdivision_factor(h)
    if m == 1:
        return word
    return "".join(ch * m for ch in word)


def material_params(word: str, lambd: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Map symbols to (D, Sigma_r, nuSigmaf).

    fissile (A / A_0): D=1.0, Sigma_r=0.5, nuSigmaf=lambda
    absorber (others): D=1.0, Sigma_r=2.0, nuSigmaf=0
    """
    symbols = np.array(list(word), dtype=object)
    is_fissile = symbols == "A"
    N = len(symbols)
    D = np.ones(N, dtype=float)
    Sigma_r = np.where(is_fissile, 0.5, 2.0).astype(float)
    nuSigmaf = np.where(is_fissile, float(lambd), 0.0).astype(float)
    return D, Sigma_r, nuSigmaf


def build_loss_matrix(D: np.ndarray, Sigma_r: np.ndarray, h: float = 1.0) -> np.ndarray:
    """
    Build symmetric tridiagonal loss operator with harmonic-mean interfaces and vacuum BC.
    """
    D = np.asarray(D, dtype=float)
    Sigma_r = np.asarray(Sigma_r, dtype=float)
    if D.ndim != 1 or Sigma_r.ndim != 1 or D.shape != Sigma_r.shape:
        raise ValueError("D and Sigma_r must be 1D arrays of equal length")
    if h <= 0:
        raise ValueError("h must be positive")

    N = D.size
    if N < 1:
        raise ValueError("empty medium")

    if N == 1:
        diag = Sigma_r[0] + 2.0 * D[0] / (h * h)
        return np.array([[diag]], dtype=float)

    D_if = 2.0 * D[:-1] * D[1:] / (D[:-1] + D[1:])
    left = np.empty(N, dtype=float)
    right = np.empty(N, dtype=float)
    left[0] = D[0]
    left[1:] = D_if
    right[:-1] = D_if
    right[-1] = D[-1]

    diag = Sigma_r + (left + right) / (h * h)
    off = -D_if / (h * h)
    L = np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)
    return L


def build_fission_matrix(nuSigmaf: np.ndarray) -> np.ndarray:
    return np.diag(np.asarray(nuSigmaf, dtype=float))


def _dominant_dense(F: np.ndarray, L: np.ndarray) -> float:
    """Return dominant real eigenvalue k of F*phi = k*L*phi using dense LAPACK."""
    vals = eig(F, L, check_finite=False, overwrite_a=False, overwrite_b=False)[0]
    vals = vals[np.isfinite(vals)]
    vals = vals[np.abs(vals.imag) < 1e-9].real
    if vals.size == 0:
        return 0.0
    return float(np.max(vals))


def _dominant_sparse(F: np.ndarray, L: np.ndarray, sigma: float = 1.0) -> float:
    """Return dominant eigenvalue using sparse shift-invert ARPACK centered at shift `sigma`."""
    F_sp = csc_matrix(F)
    L_sp = csc_matrix(L)
    vals = eigs(F_sp, M=L_sp, k=1, sigma=sigma, which="LM", return_eigenvectors=False)
    val = vals[0]
    return float(val.real)


def keff(L: np.ndarray, F: np.ndarray, sigma: float = 1.0) -> float:
    N = L.shape[0]
    # Use dense LAPACK for N <= 1000: dense is reliable and fast enough (< 1 s for N=1000).
    # The sparse path uses shift-invert ARPACK (scipy.sparse.linalg.eigs) with sigma as the
    # shift point; shift-invert solves (A - sigma*M)^{-1} M v = mu v and is designed to find
    # eigenvalues nearest to sigma.  For the generalized problem F*phi = k*L*phi the dominant
    # eigenvalue sits near 1.0 when lambda is near criticality, which is why sigma=1.0 is the
    # natural target.  However, for moderate N (400-600) the factorisation of (F-sigma*L) can
    # be ill-conditioned due to the aperiodic near-zero structure of F, causing ARPACK to
    # converge to a spurious eigenvalue.  Dense LAPACK avoids this by computing the full
    # spectrum directly, and is reserved here for truly large systems only.
    if N <= 1000:
        return _dominant_dense(F, L)
    return _dominant_sparse(F, L, sigma=sigma)


def lambda_c(
    word: str,
    bracket: tuple[float, float] = (0.3, 6.0),
    tol: float = 1e-9,
    h: float = 1.0,
) -> float:
    word_eff = refined_word_for_mesh(word, h)
    D, Sigma_r, _ = material_params(word_eff, 0.0)
    L = build_loss_matrix(D, Sigma_r, h=h)

    def f(lmbd: float) -> float:
        _, _, nuSigmaf = material_params(word_eff, lmbd)
        F = build_fission_matrix(nuSigmaf)
        return keff(L, F) - 1.0

    a, b = bracket
    fa = f(a)
    fb = f(b)
    if fa * fb > 0:
        raise ValueError(f"Root not bracketed for lambda in [{a}, {b}] (f(a)={fa}, f(b)={fb})")
    return float(brentq(f, a, b, xtol=tol, rtol=tol, maxiter=200))


@dataclass
class SmokeTestResult:
    computed_lambda_c: float
    analytic_lambda_c: float
    relative_error: float
    passed: bool


def uniform_slab_smoke_test(N: int = 50, h: float = 1.0, tol_rel: float = 0.01) -> SmokeTestResult:
    word = "A" * N
    computed = lambda_c(word, bracket=(0.3, 6.0), tol=1e-10)

    D = 1.0
    Sigma_r = 0.5
    L = N * h
    analytic = D * (math.pi / L) ** 2 + Sigma_r
    rel_err = abs(computed - analytic) / analytic
    passed = rel_err <= tol_rel
    return SmokeTestResult(
        computed_lambda_c=float(computed),
        analytic_lambda_c=float(analytic),
        relative_error=float(rel_err),
        passed=bool(passed),
    )


def run_smoke_test_or_fail() -> SmokeTestResult:
    result = uniform_slab_smoke_test()
    if not result.passed:
        raise RuntimeError(
            "Uniform slab smoke test failed: "
            f"computed={result.computed_lambda_c:.12f}, "
            f"analytic={result.analytic_lambda_c:.12f}, "
            f"relative_error={result.relative_error:.6%}"
        )
    return result


if __name__ == "__main__":
    out = run_smoke_test_or_fail()
    print(
        "Uniform-fissile slab smoke test passed:",
        f"computed λ_c={out.computed_lambda_c:.12f},",
        f"analytic λ_c={out.analytic_lambda_c:.12f},",
        f"relative error={out.relative_error:.6%}",
    )
