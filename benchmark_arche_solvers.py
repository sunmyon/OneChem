"""Repeated same-condition timing of old, unscaled, and scaled ARCHE."""

import os
from pathlib import Path
import resource
import subprocess
import time

import h5py
import numpy as np


REPO = Path("/Users/j/Documents/Codex/2026-09-10/referenced-chatgpt-conversation-this-is-an/work/arche-dev")
OLD = Path("/private/tmp/arche-old-comparison/build-release/src/apps/collapse_metal_grain/metal_collapse")
CURRENT = REPO / "results/ch-step-experiments/collapse_experiment"
ZETA = os.environ.get("ARCHE_BENCH_ZETA", "1e-17")
OUT = Path(__file__).resolve().parent / ("results/solver_speed_cr0" if ZETA == "0" else "results/solver_speed_cr1e17")
PROFILE = os.environ.get("ARCHE_BENCH_PROFILE", "0") == "1"
if PROFILE:
    OUT = Path(str(OUT) + "_kernel")
BASE_ENV = {
    "METAL_JLW21": "0", "METAL_ZETA0": ZETA, "METAL_Z_METAL": "1",
    "METAL_DT_FACTOR": "0.05", "METAL_DT_FACTOR_CHEM": "1e50",
    "METAL_DT_FACTOR_INIT": "1e-8", "METAL_N_INIT_STEPS": "10",
    "METAL_OUTPUT_STRIDE": "10", "METAL_XNH_STOP": "1e10",
    "METAL_MAX_ITER": "2000",
}


def measure(name, rep):
    folder = OUT / f"{name}_{rep}"
    folder.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update(BASE_ENV)
    env["METAL_OUTDIR"] = str(folder)
    if PROFILE:
        env["METAL_BENCH"] = "1"
    if name == "old":
        binary = OLD
        env.pop("ARCHE_SOLVER_MODE", None)
        env.pop("ARCHE_EXPERIMENT_SCALED_LINEAR", None)
    else:
        binary = CURRENT
        env.update({"ARCHE_SOLVER_MODE": "hybrid", "ARCHE_BLOCK_SOLVER": "require",
                    "ARCHE_CONTINUUM_CLOSURE": "onezone",
                    "ARCHE_EXPERIMENT_SCALED_LINEAR": "1" if name == "scaled" else "0"})
    usage_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.perf_counter()
    result = subprocess.run([str(binary)], cwd=REPO, env=env,
                            stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                            text=True)
    wall = time.perf_counter() - start
    usage_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = (usage_after.ru_utime + usage_after.ru_stime) - (usage_before.ru_utime + usage_before.ru_stime)
    path = folder / ("collapse_CR0_Z1e+0.h5" if ZETA == "0" else "collapse_CR1e-17_Z1e+0.h5")
    with h5py.File(path) as h:
        steps = int(h["step"][-1])
        density = float(h["nH"][-1])
        exit_code = int(h.attrs["exit_code"])
    kernel = None
    if PROFILE:
        bench = np.loadtxt(folder / f"bench_CR{'0' if ZETA == '0' else '1e-17'}_Z1e+0.csv",
                           delimiter=",", skiprows=1)
        kernel = float(np.sum(bench[:, 2]) * 1e-6)
    print(f"{name:8} rep={rep} wall={wall:.4f}s cpu={cpu:.4f}s "
          f"steps={steps} last_nH={density:.5e} exit={exit_code} rc={result.returncode} "
          f"kernel_wall={kernel}", flush=True)
    if result.returncode != 0 or exit_code != 0:
        print(result.stderr[-2000:])
    return wall, cpu, steps, density, exit_code, result.returncode


results = {name: [] for name in ("old", "unscaled", "scaled")}
for rep in range(1, int(os.environ.get("ARCHE_BENCH_REPS", "3")) + 1):
    order = ["old", "unscaled", "scaled"]
    order = order[rep - 1:] + order[:rep - 1]
    for name in order:
        results[name].append(measure(name, rep))

for name, values in results.items():
    a = np.asarray(values)
    print(f"SUMMARY {name:8} wall_median={np.median(a[:,0]):.4f}s "
          f"cpu_median={np.median(a[:,1]):.4f}s "
          f"cpu_range=[{np.min(a[:,1]):.4f},{np.max(a[:,1]):.4f}] "
          f"steps={int(np.median(a[:,2]))}")
