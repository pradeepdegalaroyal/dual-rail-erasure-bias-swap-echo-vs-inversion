"""Reproduce every result and figure: ``python experiments/run_all.py [--quick]``."""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from erasure_echo import figures, studies  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="small smoke run (< 1 minute)")
    ap.add_argument("--out", type=Path, default=ROOT, help="output root (results/ and figures/)")
    a = ap.parse_args()
    q = a.quick
    res, fig = a.out / "results", a.out / "figures"
    res.mkdir(parents=True, exist_ok=True)
    fig.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    seeds = [0, 1] if q else [0, 1, 2, 3, 4]
    shots, n_cal = (4000, 4000) if q else (20000, 20000)
    reps = 40 if q else 200
    n_grid = 9 if q else 25

    # H1: bias scaling
    bias = studies.bias_scaling(np.geomspace(0.05, 0.8, 8), np.array([0.1, 0.2, 0.4]))
    bias.to_csv(res / "h1_bias_scaling.csv", index=False)
    slope = {}
    for t, s in bias.groupby("hold"):
        slope[f"{t:g}"] = float(np.polyfit(np.log(s["eta"]), np.log(s["infidelity"]), 1)[0])
    max_rel_dev = float(np.max(np.abs(bias["infidelity"] / bias["analytic"] - 1)))
    with open(res / "h1_summary.json", "w") as fh:
        json.dump({"loglog_slope_vs_eta": slope, "max_rel_deviation_from_analytic": max_rel_dev},
                  fh, indent=2)

    # H2: phase maps (no detuning, isolates the non-Hermitian term)
    holds = np.geomspace(0.02, 0.8, n_grid)
    map_a = studies.map_eta_hold(np.linspace(0.05, 0.9, n_grid), holds, eps=0.02, sigma_delta=0.0)
    map_b = studies.map_eps_hold(np.geomspace(1e-4, 5e-2, n_grid), holds, eta=0.5, sigma_delta=0.0)
    map_a.to_csv(res / "h2_map_eta_hold.csv", index=False)
    map_b.to_csv(res / "h2_map_eps_hold.csv", index=False)
    figures.fig_phase_maps(map_a, map_b, fig / "fig1_phase_maps.png")

    # H3: methods comparison with shots, calibration and drift (headline holds + Pareto sweep)
    head = [0.1, 0.3, 0.6]
    pareto_holds = np.array(head) if q else np.geomspace(0.05, 0.8, 10)
    allh = np.unique(np.concatenate([pareto_holds, head]))
    meth = studies.methods_study(allh, seeds, shots, n_cal, reps, sigma_drift=0.1)
    meth.to_csv(res / "h3_methods_per_seed.csv", index=False)
    cols = ["bare", "inversion", "inversion_oracle", "echo_n2", "echo_n4", "echo_n8"]
    sub = meth[np.isclose(meth["hold"].values[:, None], head).any(axis=1)]
    summ = sub.groupby("hold")[cols].agg(["mean", "std"])
    summ.columns = [f"{c}_{s}" for c, s in summ.columns]
    summ.reset_index().to_csv(res / "h3_methods_summary.csv", index=False)
    figures.fig_dumbbell(sub, head, fig / "fig2_methods_dumbbell.png")
    ptab = studies.pareto_table(meth)
    ptab.to_csv(res / "h3_pareto_table.csv", index=False)
    figures.fig_pareto(ptab, fig / "fig3_pareto.png")

    # H3: crossover map
    nd = 5 if q else 9
    cr = studies.crossover_map(np.linspace(0.0, 0.5, nd), np.geomspace(5e-4, 3e-2, nd),
                               seeds[:3], 0.3, shots, n_cal, reps)
    cr.to_csv(res / "h3_crossover_map.csv", index=False)
    figures.fig_crossover(cr, fig / "fig4_crossover.png")

    # Ablation, pulse-count and dephasing sweeps (exact)
    studies.ablation().to_csv(res / "ablation_components.csv", index=False)
    pulses = studies.pulse_count_scaling((0, 2, 4, 8, 16, 32), (0.0, 0.005, 0.02), (0.001, 0.01))
    pulses.to_csv(res / "pulse_count_scaling.csv", index=False)
    studies.dephasing_sweep((0.0, 0.01, 0.05, 0.2), np.array([0.1, 0.3, 0.6])).to_csv(
        res / "dephasing_sweep.csv", index=False)
    figures.fig_scaling(bias, pulses, fig / "fig5_scaling.png")

    meta = {"quick": q, "seeds": seeds, "shots_per_basis": shots, "calibration_shots": n_cal,
            "reps_per_seed": reps, "runtime_s": round(time.time() - t0, 1)}
    with open(res / "run_meta.json", "w") as fh:
        json.dump(meta, fh, indent=2)
    print(f"done in {meta['runtime_s']} s; results in {res}")


if __name__ == "__main__":
    main()
