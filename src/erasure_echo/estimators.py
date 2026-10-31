"""Shot-limited logical-state estimation, with and without classical Kraus inversion."""
from __future__ import annotations

import numpy as np

from .model import PHI_MEAS, R_IN, conditional_outputs


def measurement_probs(s: np.ndarray) -> np.ndarray:
    """Outcome probabilities ``(M, 6, 3, 2)`` of logical Pauli readout (+/-) for each input."""
    rho = conditional_outputs(s)[..., 1:, 1:]
    return np.real(np.einsum("bsi,mkij,bsj->mkbs", PHI_MEAS.conj(), rho, PHI_MEAS))


def sample_bloch(s: np.ndarray, shots: int, rng: np.random.Generator) -> np.ndarray:
    """Estimate the accepted-state Bloch vector for each input from ``shots`` per basis.

    Erased shots are discarded (postselection).  Returns ``(M, 6, 3)``.
    """
    p = np.clip(measurement_probs(s), 0.0, 1.0)
    n_plus = rng.binomial(shots, p[..., 0])
    p_minus_given = np.clip(p[..., 1] / np.clip(1.0 - p[..., 0], 1e-12, None), 0.0, 1.0)
    n_minus = rng.binomial(shots - n_plus, p_minus_given)
    n_acc = n_plus + n_minus
    return np.where(n_acc > 0, (n_plus - n_minus) / np.maximum(n_acc, 1), 0.0)


def calibrate_decay(
    s_bare: np.ndarray, shots: int, rng: np.random.Generator
) -> tuple[np.ndarray, np.ndarray]:
    """Estimate ``x_a = gamma_a T`` and ``x_b = gamma_b T`` from survival of |a> and |b> inputs."""
    rho = conditional_outputs(s_bare)
    surv_a = np.real(rho[:, 4, 1, 1] + rho[:, 4, 2, 2])
    surv_b = np.real(rho[:, 5, 1, 1] + rho[:, 5, 2, 2])
    out = []
    for surv in (surv_a, surv_b):
        k = rng.binomial(shots, np.clip(surv, 0.0, 1.0))
        out.append(-np.log(np.clip(k, 0.5, None) / shots))
    return out[0], out[1]


def kraus_invert(r_hat: np.ndarray, x_a: np.ndarray, x_b: np.ndarray) -> np.ndarray:
    """Undo the no-jump map ``K = diag(exp(-x_a/2), exp(-x_b/2))`` on estimated Bloch vectors.

    ``r_hat`` has shape ``(M, 6, 3)``; the result is clipped to the Bloch ball.
    """
    ea = np.exp(x_a)[:, None]
    eb = np.exp(x_b)[:, None]
    rx, ry, rz = r_hat[..., 0], r_hat[..., 1], r_hat[..., 2]
    p_aa, p_bb = 0.5 * (1 + rz) * ea, 0.5 * (1 - rz) * eb
    coh = 0.5 * (rx - 1j * ry) * np.sqrt(ea * eb)  # rho_ab
    tr = p_aa + p_bb
    out = np.stack([2 * coh.real / tr, -2 * coh.imag / tr, (p_aa - p_bb) / tr], axis=-1)
    norm = np.linalg.norm(out, axis=-1, keepdims=True)
    return out / np.maximum(norm, 1.0)


def bloch_mse(r_est: np.ndarray) -> np.ndarray:
    """Mean squared Bloch error against the ideal inputs, per batch entry."""
    return ((r_est - R_IN[None]) ** 2).sum(axis=-1).mean(axis=-1)
