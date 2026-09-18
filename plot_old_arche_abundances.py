"""Compare density evolution of all 89 abundances with original ARCHE."""

from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


OLD = Path("/Users/j/Documents/Codex/2026-09-10/referenced-chatgpt-conversation-this-is-an/work/arche-dev/results/comparison-zsolar/dt005-nH1e16/old/collapse_CR0_Z1e+0.h5")
NEW = Path("/Users/j/Documents/Codex/2026-09-10/referenced-chatgpt-conversation-this-is-an/work/arche-dev/results/ch-step-experiments/density-n1e16-scaled/collapse_CR0_Z1e+0.h5")
OUT = Path(__file__).resolve().parent / "results" / "old_arche_comparison"
OUT.mkdir(parents=True, exist_ok=True)


def load(path):
    with h5py.File(path) as h:
        names = h["y"].attrs["species"].decode().split(",")
        return np.asarray(h["nH"]), np.asarray(h["y"]), names, dict(h.attrs)


n_old, y_old, names_old, attrs_old = load(OLD)
n_new, y_new, names_new, attrs_new = load(NEW)
assert names_old == names_new
assert y_old.shape[1] == y_new.shape[1] == 89
for key in ("Z_metal", "zeta0_cgs", "ic_T_K0", "ic_y_H2", "ic_y_HD", "ic_y_e0"):
    assert attrs_old[key] == attrs_new[key], key
assert np.all(np.diff(n_old) > 0) and np.all(np.diff(n_new) > 0)

plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9})


def draw(axis, idx, xmin=None, xmax=None):
    for n, y, label, color, line in [
        (n_old, y_old, "old ARCHE", "#d65b39", "--"),
        (n_new, y_new, "updated hybrid + scaled solve", "#196da5", "-"),
    ]:
        values = np.where(y[:, idx] > 0, y[:, idx], np.nan)
        axis.plot(n, values, label=label, color=color, lw=1.35, ls=line)
    axis.set_xscale("log")
    axis.set_yscale("log")
    if not np.any(y_old[:, idx] > 0) and not np.any(y_new[:, idx] > 0):
        axis.set_ylim(1e-40, 1e-30)
        axis.text(0.5, 0.5, "zero in both runs", transform=axis.transAxes,
                  ha="center", va="center", color="#666666")
    axis.set_title(names_old[idx], loc="left", fontweight="bold")
    axis.grid(which="major", alpha=0.22)
    axis.grid(which="minor", alpha=0.08)
    if xmin is not None:
        axis.set_xlim(xmin, xmax)
    else:
        axis.set_xlim(0.1, 1.1e16)


def make_grid(indices, output, title, columns=4, xmin=None, xmax=None):
    rows = int(np.ceil(len(indices) / columns))
    fig, axes = plt.subplots(rows, columns, figsize=(3.55 * columns, 2.55 * rows),
                             constrained_layout=True, squeeze=False)
    for axis, idx in zip(axes.flat, indices):
        draw(axis, idx, xmin, xmax)
    for axis in axes.flat[len(indices):]:
        axis.set_visible(False)
    for axis in axes[-1]:
        if axis.get_visible():
            axis.set_xlabel(r"$n_{\rm H}$ [cm$^{-3}$]")
    for axis in axes[:, 0]:
        axis.set_ylabel("abundance / H nucleus")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.035), frameon=False)
    fig.suptitle(title, y=1.075, fontsize=12)
    fig.savefig(OUT / output, dpi=180, bbox_inches="tight")
    plt.close(fig)


focus = ["H", "H2", "H+", "e-", "He", "He+", "He++", "HD", "C", "C+",
         "O", "O+", "CO", "OH", "H2O", "HCO", "CH", "CH2", "CH3", "CH4"]
make_grid([names_old.index(name) for name in focus], "selected_abundances.png",
          "Solar Z, CR=0: old ARCHE vs updated hybrid (density evolution)")
make_grid([names_old.index(name) for name in ["CH", "CH2", "CH3", "CH4"]],
          "ch_family_zoom.png", "CH-family detail, solar Z, CR=0", columns=2,
          xmin=1e7, xmax=1e10)
for start in range(0, len(names_old), 24):
    end = min(start + 24, len(names_old))
    make_grid(range(start, end), f"all_species_{start + 1:02d}_{end:02d}.png",
              f"All species {start + 1}–{end} of {len(names_old)}: old ARCHE vs updated hybrid")

print(f"old: {len(n_old)} records, final density {n_old[-1]:.6e}, exit {attrs_old['exit_message'].decode()}")
print(f"new: {len(n_new)} records, final density {n_new[-1]:.6e}, exit {attrs_new['exit_message'].decode()}")
print(OUT)
