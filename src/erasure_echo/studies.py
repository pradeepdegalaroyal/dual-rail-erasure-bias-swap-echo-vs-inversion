"""Study functions: each returns a pandas DataFrame of raw results (no plotting here)."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .estimators import bloch_mse, calibrate_decay, kraus_invert, sample_bloch
from .model import Params, channel, evaluate


def bias_scaling(etas: np.ndarray, holds: np.ndarray) -> pd.DataFrame:
    """H1: bare erasure-postselected infidelity vs asymmetry and hold time (no dephasing)."""
    rows = []
    p = Params(gamma_phi=0.0)
    for t in holds:
        inf = evaluate(channel(p, 0, t, etas))["infidelity"]
        for e, v in zip(etas, inf):
            rows.append({"eta": e, "hold": t, "infidelity": v, "analytic": (2 * e * t) ** 2 / 24})
    return pd.DataFrame(rows)


def map_eta_hold(etas: np.ndarray, holds: np.ndarray, eps: float, sigma_delta: float) -> pd.DataFrame:
    """H2 phase map: log10(infidelity_echo / infidelity_bare) over (eta, hold), echo with n = 2."""
    rows = []
    p = Params(eps=eps, sigma_delta=sigma_delta)
    for t in holds:
        bare = evaluate(channel(p, 0, t, etas))["infidelity"]
        echo = evaluate(channel(p, 2, t, etas))["infidelity"]
        for e, b, c in zip(etas, bare, echo):
            rows.append({"eta": e, "hold": t, "eps": eps, "sigma_delta": sigma_delta,
                         "infid_bare": b, "infid_echo": c, "log10_ratio": np.log10(c / b)})
    return pd.DataFrame(rows)


def map_eps_hold(epss: np.ndarray, holds: np.ndarray, eta: float, sigma_delta: float) -> pd.DataFrame:
    """H2 phase map over (pulse over-rotation, hold time) at fixed asymmetry."""
    rows = []
    for t in holds:
        for eps in epss:
            p = Params(eta=eta, eps=eps, sigma_delta=sigma_delta)
            b = evaluate(channel(p, 0, t))["infidelity"][0]
            c = evaluate(channel(p, 2, t))["infidelity"][0]
            rows.append({"eps": eps, "hold": t, "eta": eta, "sigma_delta": sigma_delta,
                         "infid_bare": b, "infid_echo": c, "log10_ratio": np.log10(c / b)})
    return pd.DataFrame(rows)


def run_methods(
    p: Params,
    hold: float,
    rng: np.random.Generator,
    shots: int,
    n_cal: int,
    reps: int,
    sigma_drift: float,
    echo_ns: tuple[int, ...] = (2, 4, 8),
    oracle: bool = True,
) -> dict[str, float]:
    """Mean shot-limited Bloch error of every method for one device draw (``reps`` repetitions).

    Loss asymmetry drifts between calibration (at ``p.eta``) and the run (``p.eta + drift``).
    """
    etas = np.clip(p.eta + sigma_drift * rng.standard_normal(reps), 0.0, 0.95)
    s_bare = channel(p, 0, hold, etas)
    r_bare = sample_bloch(s_bare, shots, rng)
    out = {"bare": bloch_mse(r_bare).mean(), "accept_bare": evaluate(s_bare)["accept"].mean()}
    s_nom = np.repeat(channel(p, 0, hold), reps, axis=0)
    xa, xb = calibrate_decay(s_nom, n_cal, rng)
    out["inversion"] = bloch_mse(kraus_invert(r_bare, xa, xb)).mean()
    if oracle:
        out["inversion_oracle"] = bloch_mse(
            kraus_invert(r_bare, (1 + etas) * hold, (1 - etas) * hold)
        ).mean()
    for n in echo_ns:
        s_echo = channel(p, n, hold, etas)
        out[f"echo_n{n}"] = bloch_mse(sample_bloch(s_echo, shots, rng)).mean()
        out[f"accept_echo_n{n}"] = evaluate(s_echo)["accept"].mean()
    return out


def methods_study(
    holds: np.ndarray, seeds: list[int], shots: int, n_cal: int, reps: int, sigma_drift: float,
    eta: float = 0.5, eps_sd: float = 0.01,
) -> pd.DataFrame:
    """Methods comparison over hold times and seeds; each seed draws its own pulse over-rotation."""
    rows = []
    for seed in seeds:
        rng = np.random.default_rng(1000 + seed)
        eps = float(rng.normal(0.0, eps_sd))
        p = Params(eta=eta, eps=eps)
        for t in holds:
            row = run_methods(p, t, rng, shots, n_cal, reps, sigma_drift)
            rows.append({"seed": seed, "hold": t, "eps": eps, "sigma_drift": sigma_drift, **row})
    return pd.DataFrame(rows)


def crossover_map(
    drifts: np.ndarray, epss: np.ndarray, seeds: list[int], hold: float, shots: int,
    n_cal: int, reps: int, eta: float = 0.5,
) -> pd.DataFrame:
    """H3 map: echo (n = 2) versus calibrated inversion over (drift, over-rotation)."""
    rows = []
    for i, d in enumerate(drifts):
        for j, e in enumerate(epss):
            inv, ech = [], []
            for seed in seeds:
                rng = np.random.default_rng(5000 + 97 * seed + 13 * i + j)
                r = run_methods(Params(eta=eta, eps=float(e)), hold, rng, shots, n_cal, reps, d,
                                echo_ns=(2,), oracle=False)
                inv.append(r["inversion"])
                ech.append(r["echo_n2"])
            rows.append({"drift": d, "eps": e, "mse_inversion": np.mean(inv), "mse_echo": np.mean(ech),
                         "log10_ratio": np.log10(np.mean(ech) / np.mean(inv))})
    return pd.DataFrame(rows)


def ablation(hold: float = 0.3) -> pd.DataFrame:
    """Deterministic ablation of echo components (conditional infidelity, exact)."""
    rows = []
    base = Params(eta=0.5, gamma_phi=0.01)
    for sd in (0.0, 1.0):
        variants = {
            "bare (postselection only)": (base.with_(sigma_delta=sd), 0),
            "echo n=2, ideal pulses": (base.with_(sigma_delta=sd, instant=True), 2),
            "echo n=2, + over-rotation 0.01": (base.with_(sigma_delta=sd, instant=True, eps=0.01), 2),
            "echo n=2, + finite pulse duration": (base.with_(sigma_delta=sd, eps=0.0), 2),
            "echo n=2, full (both)": (base.with_(sigma_delta=sd, eps=0.01), 2),
        }
        for name, (p, n) in variants.items():
            e = evaluate(channel(p, n, hold))
            rows.append({"variant": name, "sigma_delta": sd, "infidelity": e["infidelity"][0],
                         "accept": e["accept"][0]})
    return pd.DataFrame(rows)


def pulse_count_scaling(ns: tuple[int, ...], epss: tuple[float, ...], t_pis: tuple[float, ...],
                        hold: float = 0.3) -> pd.DataFrame:
    """Infidelity versus number of swap pulses for several over-rotations and pulse durations."""
    rows = []
    for eps in epss:
        for t_pi in t_pis:
            p = Params(eps=eps, t_pi=t_pi)
            for n in ns:
                if n * t_pi >= hold:
                    continue
                e = evaluate(channel(p, n, hold))
                rows.append({"eps": eps, "t_pi": t_pi, "n": n, "infidelity": e["infidelity"][0],
                             "accept": e["accept"][0]})
    return pd.DataFrame(rows)


def dephasing_sweep(gphis: tuple[float, ...], holds: np.ndarray) -> pd.DataFrame:
    """Noise sweep: where does Markovian dephasing bury the non-Hermitian bias?"""
    rows = []
    for g in gphis:
        p = Params(gamma_phi=g, eps=0.0)
        for t in holds:
            b = evaluate(channel(p, 0, t))["infidelity"][0]
            c = evaluate(channel(p.with_(instant=True), 2, t))["infidelity"][0]
            asym = evaluate(channel(p.with_(eta=0.0), 0, t))["infidelity"][0]
            rows.append({"gamma_phi": g, "hold": t, "infid_bare": b, "infid_echo_ideal": c,
                         "infid_symmetric_loss": asym, "bias_share": 1 - asym / b})
    return pd.DataFrame(rows)


def pareto_table(meth: pd.DataFrame) -> pd.DataFrame:
    """Per hold time (seed mean): rejected fraction and Bloch MSE of each method, and whether the
    method is non-dominated (no other method has both lower rejection and lower error)."""
    g = meth.groupby("hold").mean(numeric_only=True).reset_index()
    spec = [("bare", "accept_bare"), ("inversion", "accept_bare"),
            ("echo_n2", "accept_echo_n2"), ("echo_n8", "accept_echo_n8")]
    rows = []
    for _, r in g.iterrows():
        pts = [(m, 1 - r[a], r[m]) for m, a in spec]
        for m, x, y in pts:
            dom = any(x2 <= x and y2 <= y and (x2 < x or y2 < y) for m2, x2, y2 in pts if m2 != m)
            rows.append({"hold": r["hold"], "method": m, "rejected": x, "mse": y, "pareto": not dom})
    return pd.DataFrame(rows)
