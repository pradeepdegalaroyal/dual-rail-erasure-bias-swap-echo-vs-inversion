import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from erasure_echo import Params, studies  # noqa: E402


def test_methods_run_is_seeded_and_finite():
    p = Params(eps=0.005)
    a = studies.run_methods(p, 0.3, np.random.default_rng(3), 2000, 2000, 10, 0.1)
    b = studies.run_methods(p, 0.3, np.random.default_rng(3), 2000, 2000, 10, 0.1)
    assert a == b
    assert all(np.isfinite(v) for v in a.values())


def test_phase_map_sign_structure():
    df = studies.map_eps_hold(np.array([1e-4, 5e-2]), np.array([0.6]), eta=0.5, sigma_delta=0.0)
    lo, hi = df.sort_values("eps")["log10_ratio"].values
    assert lo < 0 < hi  # tiny over-rotation: echo wins; large: bare wins
