#!/usr/bin/env python3
"""Plot the one-zone parameter-validation trajectories.

Usage:
    python3 docs/plot_onezone_parameter_validation.py RAW_OUTPUT_ROOT OUTPUT_DIR
"""

from argparse import ArgumentParser
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def read_trajectory(path):
    data = np.loadtxt(path)
    return {
        "nH": 0.76 * data[:, 0],
        "temperature": data[:, 1],
        "electron": data[:, 10],
        "H2": data[:, 16],
    }


def save_figure(figure, output_dir, name):
    figure.savefig(output_dir / name, dpi=220, bbox_inches="tight")
    plt.close(figure)


def configure_axes(axes):
    for axis in np.atleast_1d(axes).flat:
        axis.grid(True, which="major", color="0.88", linewidth=0.7)
        axis.grid(True, which="minor", color="0.94", linewidth=0.4)
        axis.tick_params(direction="in", which="both", top=True, right=True)


def plot_thermal_sequences(root, output_dir):
    figure, axes = plt.subplots(1, 3, figsize=(14.0, 4.2), sharex=True)
    configure_axes(axes)

    metallicities = [(-5, r"$10^{-5}$"), (-4, r"$10^{-4}$"),
                     (-3, r"$10^{-3}$"), (-2, r"$10^{-2}$"),
                     (-1, r"$10^{-1}$"), (0, r"$1$")]
    colors = plt.cm.viridis(np.linspace(0.05, 0.95, len(metallicities)))
    base = (root / "feature_nonzero_c1_h1_d1" /
            "metals_CO_EQ_FULL_He_HD_CMB")

    for (exponent, label), color in zip(metallicities, colors):
        path = base / f"output_MT_HD_ME{exponent}_EIMP_CMBz0.dat"
        data = read_trajectory(path)
        axes[0].loglog(data["nH"], data["temperature"], color=color,
                       linewidth=1.7, label=label)

    axes[0].set_title("Metallicity")
    axes[0].set_ylabel(r"Gas temperature $T$ [K]")
    axes[0].legend(title=r"$Z/Z_\odot$", ncol=2, fontsize=8,
                   title_fontsize=9, frameon=False)

    collapse_paths = [
        ("0.1", root / "tcoll_0p1tff" / "metals_CO_EQ_FULL_He_HD_CMB" /
         "output_MT_HD_ME-3_EIMP_CMBz0_BETA_10.dat"),
        ("1", root / "tcoll_1tff" / "metals_CO_EQ_FULL_He_HD_CMB" /
         "output_MT_HD_ME-3_EIMP_CMBz0.dat"),
        ("10", root / "tcoll_10tff" / "metals_CO_EQ_FULL_He_HD_CMB" /
         "output_MT_HD_ME-3_EIMP_CMBz0_BETA_0.1.dat"),
    ]
    for label, path in collapse_paths:
        data = read_trajectory(path)
        axes[1].loglog(data["nH"], data["temperature"], linewidth=1.8,
                       label=rf"${label}\,t_{{\rm ff}}$")

    axes[1].set_title(r"Collapse time ($Z=10^{-3}Z_\odot$)")
    axes[1].legend(fontsize=8, frameon=False)

    redshifts = [0, 5, 10, 15, 20, 30]
    colors = plt.cm.plasma(np.linspace(0.05, 0.9, len(redshifts)))
    base = root / "cmb_z0_to_30" / "metals_CO_EQ_FULL_He_HD_CMB"

    for redshift, color in zip(redshifts, colors):
        path = base / f"output_MT_HD_ME-3_EIMP_CMBz{redshift}.dat"
        data = read_trajectory(path)
        axes[2].loglog(data["nH"], data["temperature"], color=color,
                       linewidth=1.7, label=rf"$z={redshift}$")

    axes[2].set_title(r"CMB redshift ($Z=10^{-3}Z_\odot$)")
    axes[2].legend(ncol=2, fontsize=8, frameon=False)

    for axis in axes:
        axis.set_xlabel(r"Hydrogen-nucleus density $n_{\rm H}$ [cm$^{-3}$]")
        axis.set_xlim(7e-4, 1.2e16)
        axis.set_ylim(2.2, 1.5e4)

    figure.tight_layout()
    save_figure(figure, output_dir, "validation_thermal_sequences.png")


def plot_chemistry_switches(root, output_dir):
    figure, axes = plt.subplots(2, 2, figsize=(11.2, 8.2))
    configure_axes(axes)

    metallicities = [(-5, r"$10^{-5}$"), (-4, r"$10^{-4}$"),
                     (-3, r"$10^{-3}$"), (-2, r"$10^{-2}$"),
                     (-1, r"$10^{-1}$"), (0, r"$1$")]
    colors = plt.cm.viridis(np.linspace(0.08, 0.92, len(metallicities)))
    co_on = (root / "feature_nonzero_c1_h1_d1" /
             "metals_CO_EQ_FULL_He_HD_CMB")
    co_off = root / "feature_nonzero_c0_h1_d1" / "metals_He_HD_CMB"
    final_on = []
    final_off = []

    for (exponent, label), color in zip(metallicities, colors):
        filename = f"output_MT_HD_ME{exponent}_EIMP_CMBz0.dat"
        on = read_trajectory(co_on / filename)
        off = read_trajectory(co_off / filename)
        axes[0, 0].loglog(on["nH"], on["temperature"], color=color,
                          linewidth=1.8)
        axes[0, 0].loglog(off["nH"], off["temperature"], color=color,
                          linewidth=1.4, linestyle="--")
        final_on.append((10.0 ** exponent, on["electron"][-1], on["H2"][-1]))
        final_off.append((10.0 ** exponent, off["electron"][-1], off["H2"][-1]))

    color_handles = [Line2D([0], [0], color=color, lw=2, label=label)
                     for (_, label), color in zip(metallicities, colors)]
    style_handles = [Line2D([0], [0], color="0.25", lw=2, label="C/O on"),
                     Line2D([0], [0], color="0.25", lw=1.5, linestyle="--",
                            label="C/O off")]
    legend1 = axes[0, 0].legend(handles=color_handles, title=r"$Z/Z_\odot$",
                                fontsize=8, frameon=False, loc="lower left")
    axes[0, 0].add_artist(legend1)
    axes[0, 0].legend(handles=style_handles, fontsize=8, frameon=False,
                      loc="upper right")
    axes[0, 0].set_title("C/O chemistry switch")
    axes[0, 0].set_ylabel(r"Gas temperature $T$ [K]")

    feature_styles = [
        (0, 0, "He off, HD off", "#4c78a8", "-"),
        (0, 1, "He off, HD on", "#f58518", "--"),
        (1, 0, "He on, HD off", "#54a24b", "-."),
        (1, 1, "He on, HD on", "#b279a2", ":"),
    ]
    for helium, hd, label, color, linestyle in feature_styles:
        directory = (root / f"feature_zero_c0_h{helium}_d{hd}" /
                     ("metals" + ("_He" if helium else "") +
                      ("_HD" if hd else "") + "_CMB"))
        filename = ("output_MT" + ("_HD" if hd else "") +
                    "_ME-inf_EIMP_CMBz0.dat")
        data = read_trajectory(directory / filename)
        axes[0, 1].loglog(data["nH"], data["temperature"], label=label,
                          color=color, linestyle=linestyle, linewidth=1.8)

    axes[0, 1].set_title("He and HD switches ($Z=0$, C/O off)")
    axes[0, 1].legend(fontsize=8, frameon=False)

    final_on = np.asarray(final_on)
    final_off = np.asarray(final_off)
    axes[1, 0].loglog(final_on[:, 0], final_on[:, 1], "o-", label="C/O on")
    axes[1, 0].loglog(final_off[:, 0], final_off[:, 1], "s--", label="C/O off")
    axes[1, 0].set_xlabel(r"Metallicity $Z/Z_\odot$")
    axes[1, 0].set_ylabel(r"Final electron abundance $x_e$")
    axes[1, 0].legend(fontsize=8, frameon=False)

    axes[1, 1].semilogx(final_on[:, 0], final_on[:, 2], "o-", label="C/O on")
    axes[1, 1].semilogx(final_off[:, 0], final_off[:, 2], "s--", label="C/O off")
    axes[1, 1].set_xlabel(r"Metallicity $Z/Z_\odot$")
    axes[1, 1].set_ylabel(r"Final molecular abundance $x_{\rm H_2}$")
    axes[1, 1].ticklabel_format(axis="y", style="plain", useOffset=False)
    axes[1, 1].legend(fontsize=8, frameon=False)

    for axis in axes[0, :]:
        axis.set_xlabel(r"Hydrogen-nucleus density $n_{\rm H}$ [cm$^{-3}$]")
        axis.set_xlim(7e-4, 1.2e16)
        axis.set_ylim(2.2, 1.5e4)

    figure.tight_layout()
    save_figure(figure, output_dir, "validation_chemistry_switches.png")


def plot_lw_sequences(root, output_dir):
    figure, axes = plt.subplots(2, 2, figsize=(11.5, 8.2), sharex=True)
    configure_axes(axes)

    j21_values = [0, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000, 10000]
    labels = ["0", r"$10^{-3}$", r"$10^{-2}$", r"$10^{-1}$", "1",
              "10", r"$10^2$", r"$10^3$", r"$10^4$"]
    colors = plt.cm.turbo(np.linspace(0.02, 0.94, len(j21_values)))
    cases = [
        ("j21_z1em3", "metals_CO_EQ_FULL_He_HD_LW_CMB", "-3",
         r"$Z=10^{-3}Z_\odot$"),
        ("j21_z0", "metals_He_HD_LW_CMB", "-inf", r"$Z=0$ (C/O off)"),
    ]

    for column, (suite, directory, metallicity, title) in enumerate(cases):
        base = root / suite / directory
        for j21, label, color in zip(j21_values, labels, colors):
            jtext = f"{j21:g}"
            path = base / f"output_MT_HD_ME{metallicity}_EIMP_B_{jtext}_CMBz0.dat"
            data = read_trajectory(path)
            axes[0, column].loglog(data["nH"], data["temperature"],
                                   color=color, linewidth=1.55, label=label)
            axes[1, column].loglog(data["nH"], data["H2"], color=color,
                                   linewidth=1.55)

        axes[0, column].set_title(title)
        axes[1, column].set_xlabel(
            r"Hydrogen-nucleus density $n_{\rm H}$ [cm$^{-3}$]")
        axes[0, column].set_xlim(7e-4, 1.2e16)
        axes[1, column].set_xlim(7e-4, 1.2e16)
        axes[0, column].set_ylim(20, 1.5e4)
        axes[1, column].set_ylim(1e-12, 0.8)

    axes[0, 0].set_ylabel(r"Gas temperature $T$ [K]")
    axes[1, 0].set_ylabel(r"Molecular abundance $x_{\rm H_2}$")
    axes[0, 1].legend(title=r"$J_{21}$", ncol=3, fontsize=8,
                      title_fontsize=9, frameon=False)

    figure.tight_layout()
    save_figure(figure, output_dir, "validation_lw_sequences.png")


def plot_zero_metal_failure(root, output_dir):
    figure, axes = plt.subplots(1, 3, figsize=(13.5, 4.0), sharex=True)
    configure_axes(axes)

    cases = [
        ("C/O off: evolving", root / "feature_zero_c0_h1_d1" /
         "metals_He_HD_CMB" / "output_MT_HD_ME-inf_EIMP_CMBz0.dat",
         "#4c78a8", "-"),
        ("C/O on: failed", root / "feature_zero_c1_h1_d1" /
         "metals_CO_EQ_FULL_He_HD_CMB" /
         "output_MT_HD_ME-inf_EIMP_CMBz0.dat", "#e45756", "--"),
    ]

    for label, path, color, linestyle in cases:
        data = read_trajectory(path)
        axes[0].loglog(data["nH"], data["temperature"], label=label,
                       color=color, linestyle=linestyle, linewidth=2.0)
        axes[1].loglog(data["nH"], data["electron"], color=color,
                       linestyle=linestyle, linewidth=2.0)
        axes[2].loglog(data["nH"], data["H2"], color=color,
                       linestyle=linestyle, linewidth=2.0)

    axes[0].set_ylabel(r"Gas temperature $T$ [K]")
    axes[1].set_ylabel(r"Electron abundance $x_e$")
    axes[2].set_ylabel(r"Molecular abundance $x_{\rm H_2}$")
    axes[0].legend(fontsize=9, frameon=False)

    for axis in axes:
        axis.set_xlabel(r"Hydrogen-nucleus density $n_{\rm H}$ [cm$^{-3}$]")
        axis.set_xlim(7e-4, 1.2e16)

    axes[0].set_ylim(70, 1.3e4)
    axes[1].set_ylim(1e-13, 1e-3)
    axes[2].set_ylim(1e-7, 0.8)
    figure.suptitle("Exact-zero-metal validation (He and HD on)", y=1.02)
    figure.tight_layout()
    save_figure(figure, output_dir, "validation_zero_metal_failure.png")


def main():
    parser = ArgumentParser()
    parser.add_argument("raw_output_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    arguments = parser.parse_args()

    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    plot_thermal_sequences(arguments.raw_output_root, arguments.output_dir)
    plot_chemistry_switches(arguments.raw_output_root, arguments.output_dir)
    plot_lw_sequences(arguments.raw_output_root, arguments.output_dir)
    plot_zero_metal_failure(arguments.raw_output_root, arguments.output_dir)


if __name__ == "__main__":
    main()
