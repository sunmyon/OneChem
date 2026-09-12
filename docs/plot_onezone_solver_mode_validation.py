#!/usr/bin/env python3
"""Plot the 2026-09-12 one-zone solver-mode validation."""

from argparse import ArgumentParser
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


METALLICITIES = [(-5, r"$10^{-5}$"), (-4, r"$10^{-4}$"),
                 (-3, r"$10^{-3}$"), (-2, r"$10^{-2}$"),
                 (-1, r"$10^{-1}$"), (0, r"$1$")]


def read_trajectory(path):
    data = np.loadtxt(path, usecols=(0, 1, 10, 16))
    return 0.76 * data[:, 0], data[:, 1], data[:, 2], data[:, 3]


def trajectory(root, mode, exponent, tag):
    matches = list((root / mode).glob(
        f"**/output_MT_HD_ME{exponent}_{tag}_CMBz0.dat"))
    if len(matches) != 1:
        raise RuntimeError(f"expected one trajectory, found {matches}")
    return read_trajectory(matches[0])


def configure(axis):
    axis.grid(True, which="major", color="0.87", linewidth=0.7)
    axis.grid(True, which="minor", color="0.94", linewidth=0.4)
    axis.tick_params(direction="in", which="both", top=True, right=True)


def plot_common_range(short_root, output_dir):
    modes = [
        ("explicit", "EXP", "explicit", "#e45756", "--"),
        ("chemistry_implicit", "CIMP", "chemistry implicit", "#72b7b2", "-."),
        ("energy_implicit", "EIMP", "energy implicit", "#4c78a8", "-"),
        ("hybrid", "HYB", "hybrid", "#f2cf5b", ":"),
    ]
    figure, axes = plt.subplots(2, 3, figsize=(13.0, 8.0), sharex=True,
                                sharey=True)

    for axis, (exponent, label) in zip(axes.flat, METALLICITIES):
        configure(axis)
        for mode, tag, mode_label, color, linestyle in modes:
            density, temperature, _, _ = trajectory(
                short_root, mode, exponent, tag)
            axis.loglog(density, temperature, color=color,
                        linestyle=linestyle, linewidth=1.8, label=mode_label)
        axis.set_title(rf"$Z/Z_\odot={label[1:-1]}$")
        axis.set_xlim(7e-4, 1.1e10)
        axis.set_ylim(2.2, 1.5e4)

    for axis in axes[-1, :]:
        axis.set_xlabel(r"Hydrogen-nucleus density $n_{\rm H}$ [cm$^{-3}$]")
    for axis in axes[:, 0]:
        axis.set_ylabel(r"Gas temperature $T$ [K]")

    axes[0, 0].legend(fontsize=8, frameon=False)
    figure.suptitle(r"Solver-mode comparison to $n_{\rm H}=10^{10}$ cm$^{-3}$")
    figure.tight_layout()
    figure.savefig(output_dir / "solver_modes_temperature_n1e10.png", dpi=220,
                   bbox_inches="tight")
    plt.close(figure)


def plot_full_implicit(full_root, output_dir):
    modes = [
        ("energy_implicit", "EIMP", "energy implicit", "#4c78a8", "-"),
        ("hybrid", "HYB", "hybrid", "#f2cf5b", "--"),
        ("chemistry_implicit", "CIMP", "chemistry implicit", "#72b7b2", "-."),
    ]
    figure, axes = plt.subplots(2, 3, figsize=(13.0, 8.0), sharex=True,
                                sharey=True)

    for axis, (exponent, label) in zip(axes.flat, METALLICITIES):
        configure(axis)
        for mode, tag, mode_label, color, linestyle in modes:
            density, temperature, _, _ = trajectory(
                full_root, mode, exponent, tag)
            valid = np.isfinite(temperature) & (temperature > 0.0)
            axis.loglog(density[valid], temperature[valid], color=color,
                        linestyle=linestyle, linewidth=1.8, label=mode_label)
        axis.set_title(rf"$Z/Z_\odot={label[1:-1]}$")
        axis.set_xlim(7e-4, 1.2e16)
        axis.set_ylim(2.2, 1.5e4)

        if exponent == 0:
            axis.text(0.97, 0.05, "CIMP: negative T / NaN",
                      color="#2a9d8f", ha="right", va="bottom",
                      transform=axis.transAxes, fontsize=8)

    for axis in axes[-1, :]:
        axis.set_xlabel(r"Hydrogen-nucleus density $n_{\rm H}$ [cm$^{-3}$]")
    for axis in axes[:, 0]:
        axis.set_ylabel(r"Gas temperature $T$ [K]")

    axes[0, 0].legend(fontsize=8, frameon=False)
    figure.suptitle("Implicit-mode full-collapse comparison")
    figure.tight_layout()
    figure.savefig(output_dir / "solver_modes_temperature_full.png", dpi=220,
                   bbox_inches="tight")
    plt.close(figure)


def main():
    parser = ArgumentParser()
    parser.add_argument("validation_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    plot_common_range(args.validation_root / "output_n1e10_serial",
                      args.output_dir)
    plot_full_implicit(args.validation_root / "plot_full", args.output_dir)


if __name__ == "__main__":
    main()
