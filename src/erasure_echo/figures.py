"""Figure generation (viridis palette, 200 dpi)."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DPI = 200
V = plt.colormaps["viridis"]


def _grid(df: pd.DataFrame, x: str, y: str, z: str):
    piv = df.pivot_table(index=y, columns=x, values=z)
    return piv.columns.values, piv.index.values, piv.values


def _contour(ax, df, x, y, z, xlabel, ylabel, title, logx=False, logy=False):
    xs, ys, zs = _grid(df, x, y, z)
    lim = max(abs(np.nanmin(zs)), abs(np.nanmax(zs)), 0.5)
    levels = np.linspace(-lim, lim, 15)
    cs = ax.contourf(xs, ys, zs, levels=levels, cmap="viridis")
    if np.nanmin(zs) < 0 < np.nanmax(zs):
        ax.contour(xs, ys, zs, levels=[0.0], colors="white", linewidths=2.2)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    if logx:
        ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")
    return cs


def fig_phase_maps(map_a: pd.DataFrame, map_b: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
    cs = _contour(axes[0], map_a, "eta", "hold", "log10_ratio",
                  "loss asymmetry eta (dimensionless)", "hold time (mean lifetimes)",
                  "(a) pulse over-rotation 0.02, no detuning", logy=True)
    _contour(axes[1], map_b, "eps", "hold", "log10_ratio",
             "swap-pulse over-rotation (fraction of pi)", "hold time (mean lifetimes)",
             "(b) loss asymmetry eta = 0.5, no detuning", logx=True, logy=True)
    fig.colorbar(cs, ax=axes, label="log10(infidelity echo / infidelity bare)  (<0: echo wins)")
    fig.suptitle("Where a two-pulse swap echo beats bare erasure postselection", fontsize=11)
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def fig_dumbbell(df: pd.DataFrame, holds: list[float], path: Path) -> None:
    fig, axes = plt.subplots(1, len(holds), figsize=(4.2 * len(holds), 4), sharey=True)
    axes = np.atleast_1d(axes)
    cols = {"inversion": V(0.15), "echo_n2": V(0.8)}
    for ax, t in zip(axes, holds):
        sub = df[np.isclose(df["hold"], t)].sort_values("seed")
        y = np.arange(len(sub))
        for k, (_, r) in enumerate(sub.iterrows()):
            for m, c in cols.items():
                ax.plot([r["bare"], r[m]], [k + (0.12 if m == "echo_n2" else -0.12)] * 2, color=c, lw=1.5)
                ax.scatter(r[m], k + (0.12 if m == "echo_n2" else -0.12), color=c, zorder=3, s=28)
            ax.scatter(r["bare"], k, color="black", zorder=3, s=22)
        ax.set_xscale("log")
        ax.set_yticks(y)
        ax.set_yticklabels([f"seed {s}" for s in sub["seed"]])
        ax.set_xlabel("Bloch-vector mean squared error")
        ax.set_title(f"hold = {t:g} lifetimes", fontsize=10)
    axes[0].scatter([], [], color="black", s=22, label="bare postselection")
    axes[0].scatter([], [], color=cols["inversion"], s=28, label="calibrated inversion")
    axes[0].scatter([], [], color=cols["echo_n2"], s=28, label="swap echo, n = 2")
    fig.legend(*axes[0].get_legend_handles_labels(), fontsize=8, loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, -0.08))
    fig.suptitle("Per-device change from bare postselection (black) to each mitigation", fontsize=11)
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def fig_pareto(tab: pd.DataFrame, path: Path) -> None:
    """Scatter of rejected fraction (cost) vs Bloch error (quality); rings mark methods that are
    non-dominated at their hold time (each hold time is a separate decision)."""
    labels = {"bare": ("bare postselection", "o"), "inversion": ("calibrated inversion", "s"),
              "echo_n2": ("swap echo n = 2", "^"), "echo_n8": ("swap echo n = 8", "v")}
    fig, ax = plt.subplots(figsize=(6.6, 4.7))
    for k, (m, (label, mk)) in enumerate(labels.items()):
        s = tab[tab["method"] == m]
        ax.scatter(s["rejected"], s["mse"], marker=mk, color=V(0.1 + 0.27 * k), s=36, label=label, alpha=0.9)
    front = tab[tab["pareto"]]
    ax.scatter(front["rejected"], front["mse"], facecolors="none", edgecolors="black", s=100,
               label="non-dominated at its hold time")
    ax.set_yscale("log")
    ax.set_xlabel("rejected fraction of shots (erasure flags plus losses)")
    ax.set_ylabel("Bloch-vector mean squared error")
    ax.set_title("Cost-quality trade-off across hold times (seed mean)", fontsize=10)
    ax.legend(fontsize=8)
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def fig_crossover(cr: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    cs = _contour(ax, cr, "drift", "eps", "log10_ratio",
                  "calibration drift of loss asymmetry (std of eta)",
                  "swap-pulse over-rotation (fraction of pi)",
                  "Swap echo versus calibrated inversion (hold = 0.3 lifetimes)", logy=True)
    fig.colorbar(cs, ax=ax, label="log10(MSE echo / MSE inversion)  (<0: echo wins)")
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def fig_scaling(bias: pd.DataFrame, pulses: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    holds = sorted(bias["hold"].unique())
    for k, t in enumerate(holds):
        s = bias[bias["hold"] == t]
        c = V(0.1 + 0.8 * k / max(len(holds) - 1, 1))
        axes[0].plot(s["eta"], s["infidelity"], "o", color=c, label=f"hold {t:g}")
        axes[0].plot(s["eta"], s["analytic"], "-", color=c, lw=1)
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_xticks([0.05, 0.1, 0.2, 0.4, 0.8])
    axes[0].set_xticklabels(["0.05", "0.1", "0.2", "0.4", "0.8"])
    axes[0].minorticks_off()
    axes[0].set_xlabel("loss asymmetry eta (dimensionless)")
    axes[0].set_ylabel("conditional infidelity after erasure postselection")
    axes[0].set_title("(a) non-Hermitian bias: simulation (points) vs (2 eta T)^2 / 24 (lines)", fontsize=9)
    axes[0].legend(fontsize=8, title="mean lifetimes")
    combos = sorted(set(zip(pulses["eps"], pulses["t_pi"])))
    for k, (eps, t_pi) in enumerate(combos):
        s = pulses[(pulses["eps"] == eps) & (pulses["t_pi"] == t_pi)].sort_values("n")
        axes[1].plot(s["n"], s["infidelity"], "o-", color=V(0.05 + 0.9 * k / max(len(combos) - 1, 1)),
                     label=f"eps {eps:g}, t_pi {t_pi:g}")
    axes[1].set_xscale("log", base=2)
    axes[1].set_yscale("log")
    axes[1].set_xlabel("number of swap pulses n in the hold")
    axes[1].set_ylabel("conditional infidelity")
    axes[1].set_title("(b) pulse-count scaling at hold = 0.3, eta = 0.5", fontsize=9)
    axes[1].legend(fontsize=7)
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
