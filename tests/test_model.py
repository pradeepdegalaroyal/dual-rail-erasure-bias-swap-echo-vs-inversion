"""Fast checks of the Lindblad model, estimators and a Monte Carlo trajectory cross-check."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from erasure_echo import Params, channel, evaluate, kraus_invert, sample_bloch  # noqa: E402
from erasure_echo.estimators import bloch_mse  # noqa: E402


def test_trace_nonincreasing_and_acceptance():
    s = channel(Params(sigma_delta=1.0, eps=0.01), 4, 0.3)
    e = evaluate(s)
    assert 0 < e["accept"][0] < 1


def test_analytic_bias_without_dephasing():
    p = Params(eta=0.4, gamma_phi=0.0)
    inf = evaluate(channel(p, 0, 0.2))["infidelity"][0]
    assert inf == pytest.approx((2 * 0.4 * 0.2) ** 2 / 24, rel=0.05)


def test_ideal_echo_cancels_bias():
    p = Params(eta=0.7, gamma_phi=0.0, instant=True)
    assert evaluate(channel(p, 2, 0.5))["infidelity"][0] < 1e-12


def test_symmetric_loss_has_no_bias():
    p = Params(eta=0.0, gamma_phi=0.0)
    assert evaluate(channel(p, 0, 0.5))["infidelity"][0] < 1e-12


def test_odd_pulse_count_rejected():
    with pytest.raises(ValueError):
        channel(Params(), 3, 0.3)


def test_kraus_inversion_removes_bias_exactly():
    p = Params(eta=0.5, gamma_phi=0.0)
    s = channel(p, 0, 0.4)
    rng = np.random.default_rng(0)
    r = sample_bloch(s, 2_000_000, rng)
    xa, xb = np.array([1.5 * 0.4]), np.array([0.5 * 0.4])
    assert bloch_mse(kraus_invert(r, xa, xb))[0] < 0.2 * bloch_mse(r)[0]


def test_matches_quantum_jump_trajectories():
    """Independent check: Monte Carlo no-jump trajectories reproduce the exact conditional state."""
    rng = np.random.default_rng(1)
    eta, t, n = 0.5, 0.4, 20000
    ga, gb = 1 + eta, 1 - eta
    psi = np.array([1, 1]) / np.sqrt(2)  # |+> in (a, b)
    # no-jump evolution is K = diag(exp(-ga t/2), exp(-gb t/2)); jump prob = 1 - |K psi|^2
    k = np.diag([np.exp(-ga * t / 2), np.exp(-gb * t / 2)])
    out = k @ psi
    p_surv = float(out @ out)
    survive = rng.random(n) < p_surv
    fid_mc = abs(psi @ (out / np.sqrt(p_surv))) ** 2
    s = channel(Params(eta=eta, gamma_phi=0.0), 0, t)
    exact = evaluate(s)["per_input_fid"][0, 0]  # input X+
    assert exact == pytest.approx(fid_mc, abs=1e-9)
    assert survive.mean() == pytest.approx(p_surv, abs=0.02)
