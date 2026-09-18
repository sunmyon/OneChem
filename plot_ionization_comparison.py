"""Plot matched ionizing-background comparisons of old and updated ARCHE."""

from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent / "results" / "ionization_comparison"
RUNS = [
    ("CR: zeta=1e-17 s^-1", "cr_old", "cr_updated", "collapse_CR1e-17_Z1e+0.h5"),
    ("X-ray: zeta_X=1e-20 s^-1, E_X=300 eV", "xray_old_1e20",
     "xray_updated_1e20", "collapse_CR0_Z1e+0_ZX1e-20.h5"),
]


def load(directory, filename):
    with h5py.File(ROOT / directory / filename) as h:
        names = h["y"].attrs["species"].decode().split(",")
        return {
            "n": h["nH"][:], "y": h["y"][:], "names": names,
            "T": h["T_K"][:], "charge": h["charge_imbal"][:],
            "attrs": dict(h.attrs),
        }


fig, axes = plt.subplots(2, 2, figsize=(11, 7.2), sharex=True,
                         constrained_layout=True)
fig_t, axes_t = plt.subplots(1, 2, figsize=(10.5, 3.6), sharex=True,
                             constrained_layout=True)
for row, (title, old_dir, new_dir, filename) in enumerate(RUNS):
    old = load(old_dir, filename)
    new = load(new_dir, filename)
    assert old["names"] == new["names"]
    for col, species in enumerate(("H+", "e-")):
        axis = axes[row, col]
        index = old["names"].index(species)
        for run, label, color, style in (
            (old, "old ARCHE", "#d45739", "--"),
            (new, "updated hybrid + scaled solve", "#156fa2", "-"),
        ):
            yy = np.where(run["y"][:, index] > 0, run["y"][:, index], np.nan)
            axis.plot(run["n"], yy, label=label, c=color, ls=style, lw=1.6)
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlim(1e4, 1.1e10)
        axis.set_title(f"{title} — {species}", loc="left", fontsize=10)
        axis.grid(alpha=0.18, which="both")
        if col == 0:
            axis.set_ylabel("abundance / H nucleus")
        if row == 1:
            axis.set_xlabel(r"$n_{\rm H}$ [cm$^{-3}$]")
    axis_t = axes_t[row]
    for run, label, color, style in (
        (old, "old ARCHE", "#d45739", "--"),
        (new, "updated hybrid + scaled solve", "#156fa2", "-"),
    ):
        axis_t.plot(run["n"], run["T"], label=label, c=color, ls=style, lw=1.5)
    axis_t.set_xscale("log")
    axis_t.set_xlim(1e4, 1.1e10)
    axis_t.set_title(title, loc="left", fontsize=10)
    axis_t.set_xlabel(r"$n_{\rm H}$ [cm$^{-3}$]")
    axis_t.set_ylabel("temperature [K]")
    axis_t.grid(alpha=0.18, which="both")

handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False,
           bbox_to_anchor=(0.5, 1.055))
fig.suptitle("Ionizing backgrounds: old ARCHE vs updated hybrid", y=1.09)
fig.savefig(ROOT / "ion_abundances.png", dpi=190, bbox_inches="tight")
fig_t.legend(handles, labels, loc="upper center", ncol=2, frameon=False,
             bbox_to_anchor=(0.5, 1.13))
fig_t.suptitle("Temperature under ionizing backgrounds", y=1.19)
fig_t.savefig(ROOT / "temperatures.png", dpi=190, bbox_inches="tight")
print(ROOT / "ion_abundances.png")
