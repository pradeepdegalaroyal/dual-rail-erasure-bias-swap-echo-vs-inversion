"""Lindblad model of a dual-rail erasure qubit with asymmetric photon loss.

Hilbert space (3 levels): ``|0>`` = both modes empty (erased), ``|a>``/``|b>`` = single
excitation in mode a / b (logical computational basis).  Time is measured in units of the
mean excitation lifetime ``1 / gamma_bar``, so ``gamma_bar * T`` is the hold time in
lifetimes.  The loss rates are ``gamma_a = 1 + eta`` and ``gamma_b = 1 - eta``.

Every sequence is returned as a 9 x 9 row-major superoperator acting on ``vec(rho)``.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy.linalg import expm

DIM = 3
_I3 = np.eye(DIM)


def _ket(i: int) -> np.ndarray:
    v = np.zeros(DIM, dtype=complex)
    v[i] = 1.0
    return v


def _op(i: int, j: int) -> np.ndarray:
    return np.outer(_ket(i), _ket(j).conj())


def _ham_super(h: np.ndarray) -> np.ndarray:
    """Superoperator of ``rho -> -i [h, rho]`` (row-major vectorisation)."""
    return -1j * (np.kron(h, _I3) - np.kron(_I3, h.T))


def _diss_super(l_op: np.ndarray) -> np.ndarray:
    """Superoperator of the Lindblad dissipator D[l_op]."""
    ldl = l_op.conj().T @ l_op
    return np.kron(l_op, l_op.conj()) - 0.5 * np.kron(ldl, _I3) - 0.5 * np.kron(_I3, ldl.T)


# Basis superoperators; the Liouvillian is linear in the rates and Hamiltonian parameters.
_L_LOSS_A = _diss_super(_op(0, 1))
_L_LOSS_B = _diss_super(_op(0, 2))
_L_PHI = _diss_super(_op(1, 1)) + _diss_super(_op(2, 2))
_L_DELTA = _ham_super(0.5 * (_op(1, 1) - _op(2, 2)))
_L_COUP = _ham_super(0.5 * (_op(1, 2) + _op(2, 1)))

# Logical input states (Pauli eigenstates) embedded in the 3-level space, order X+ X- Y+ Y- Z+ Z-.
_LOGICAL = np.array(
    [[1, 1], [1, -1], [1, 1j], [1, -1j], [1, 0], [0, 1]], dtype=complex
)
_LOGICAL[:4] /= np.sqrt(2)
PSI_IN = np.concatenate([np.zeros((6, 1), dtype=complex), _LOGICAL], axis=1)
R_IN = np.array(
    [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], dtype=float
)
# Measurement eigenstates per Pauli axis (X, Y, Z) and sign (+, -): shape (3, 2, 2).
PHI_MEAS = np.array(
    [
        [[1, 1], [1, -1]],
        [[1, 1j], [1, -1j]],
        [[1, 0], [0, 1]],
    ],
    dtype=complex,
)
PHI_MEAS[:2] /= np.sqrt(2)
RHO_IN_VEC = np.array([np.outer(p, p.conj()).reshape(-1) for p in PSI_IN])  # (6, 9)


@dataclass(frozen=True)
class Params:
    """Device and control parameters (all rates in units of gamma_bar)."""

    eta: float = 0.5  # loss-rate asymmetry, gamma_{a,b} = 1 +/- eta
    gamma_phi: float = 0.01  # differential dephasing rate
    sigma_delta: float = 0.0  # std of quasi-static mode detuning
    eps: float = 0.0  # fractional over-rotation of every swap pulse
    t_pi: float = 0.003  # swap-pulse duration
    instant: bool = False  # True: ideal instantaneous (lossless) swap pulses

    def with_(self, **kw) -> "Params":
        return replace(self, **kw)


def _liouvillian(eta: np.ndarray, delta: np.ndarray, gamma_phi: float, coup: float) -> np.ndarray:
    eta = np.asarray(eta, dtype=float)
    delta = np.asarray(delta, dtype=float)
    return (
        (1.0 + eta)[:, None, None] * _L_LOSS_A
        + (1.0 - eta)[:, None, None] * _L_LOSS_B
        + gamma_phi * _L_PHI
        + delta[:, None, None] * _L_DELTA
        + coup * _L_COUP
    )


def _swap_unitary_super(theta: float) -> np.ndarray:
    """Instantaneous logical X rotation by angle ``theta`` (pi = perfect swap)."""
    u = np.eye(DIM, dtype=complex)
    c, s = np.cos(theta / 2), -1j * np.sin(theta / 2)
    u[1:, 1:] = np.array([[c, s], [s, c]])
    return np.kron(u, u.conj())


def _delta_nodes(sigma: float, n_nodes: int = 7) -> tuple[np.ndarray, np.ndarray]:
    if sigma <= 0:
        return np.zeros(1), np.ones(1)
    x, w = np.polynomial.hermite_e.hermegauss(n_nodes)
    return sigma * x, w / np.sqrt(2 * np.pi)


def channel(
    p: Params, n_pulses: int, t_total: float, etas: np.ndarray | None = None
) -> np.ndarray:
    """Superoperators of a hold of total duration ``t_total`` with ``n_pulses`` swap pulses.

    ``etas`` (shape ``(M,)``) overrides ``p.eta`` per batch entry (used for drift).  The
    quasi-static detuning is averaged exactly over a Gauss-Hermite rule.  Returns ``(M, 9, 9)``.
    The pulse train is ``tau/2, pi, tau, pi, ..., pi, tau/2`` (CPMG spacing).
    """
    if n_pulses % 2:
        raise ValueError("n_pulses must be even so the net logical operation is the identity")
    etas = np.array([p.eta]) if etas is None else np.atleast_1d(np.asarray(etas, dtype=float))
    m = etas.size
    d_nodes, d_w = _delta_nodes(p.sigma_delta)
    out = np.zeros((m, DIM**2, DIM**2), dtype=complex)
    theta = np.pi * (1.0 + p.eps)
    for dj, wj in zip(d_nodes, d_w):
        dl = np.full(m, dj)
        if n_pulses == 0:
            s = expm(_liouvillian(etas, dl, p.gamma_phi, 0.0) * t_total)
        else:
            tau = (t_total - (0.0 if p.instant else n_pulses * p.t_pi)) / n_pulses
            if tau <= 0:
                raise ValueError("pulses do not fit in the hold time")
            lh = expm(_liouvillian(etas, dl, p.gamma_phi, 0.0) * (tau / 2))
            lf = expm(_liouvillian(etas, dl, p.gamma_phi, 0.0) * tau)
            if p.instant:
                pulse = np.broadcast_to(_swap_unitary_super(theta), (m, 9, 9))
            else:
                coup = theta / p.t_pi
                pulse = expm(_liouvillian(etas, dl, p.gamma_phi, coup) * p.t_pi)
            s = lh
            for k in range(n_pulses):
                s = pulse @ s
                s = (lf if k < n_pulses - 1 else lh) @ s
        out += wj * s
    return out


def conditional_outputs(s: np.ndarray) -> np.ndarray:
    """Unnormalised 3 x 3 output states for the six Pauli inputs, shape ``(M, 6, 3, 3)``."""
    out = np.einsum("mij,kj->mki", s, RHO_IN_VEC)
    return out.reshape(s.shape[0], 6, DIM, DIM)


def evaluate(s: np.ndarray) -> dict[str, np.ndarray]:
    """Erasure-postselected figures of merit for a batch of superoperators.

    ``infidelity`` is the mean over the six Pauli inputs of ``1 - <psi|rho_cond|psi>`` with
    ``rho_cond`` the accepted (non-erased) state renormalised; ``accept`` is the mean acceptance
    probability.
    """
    rho = conditional_outputs(s)
    acc = np.real(rho[..., 1, 1] + rho[..., 2, 2])  # (M, 6)
    psi = PSI_IN[:, 1:]
    num = np.real(np.einsum("ki,mkij,kj->mk", psi.conj(), rho[..., 1:, 1:], psi))
    fid = num / acc
    return {"infidelity": 1.0 - fid.mean(axis=1), "accept": acc.mean(axis=1), "per_input_fid": fid}
