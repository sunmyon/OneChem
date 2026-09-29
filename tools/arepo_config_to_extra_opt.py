#!/usr/bin/env python3
"""Extract the ONEZONE-relevant compile options from an AREPO Config.sh.

Problem this solves: ONEZONE and AREPO's chemistry_ONECHEM module are meant
to share chemistry_solver.cc byte-for-byte, and their behavior is controlled
by preprocessor macros (#ifdef ...).  Historically the ONEZONE side of a
comparison build was assembled by hand (an ad-hoc EXTRA_OPT string typed on
the `make` command line), which silently drifted out of sync with whatever
Config.sh actually produced a given AREPO run -- exactly the kind of gap
that makes replaying a captured chemistry failure solve a subtly different
problem than the real run did.

This script removes the hand-typing step: it derives the whitelist of
"chemistry-relevant" macro names directly from ONEZONE's own source (every
identifier that appears in a #ifdef/#ifndef/#if defined(...) anywhere under
src/chemistry, src/radiation, src/standalone -- deliberately excluding
legacy/), then intersects that whitelist with whatever options are actually
active in the given AREPO Config.sh. This is self-maintaining: if a new
#ifdef FOO appears in ONEZONE's source and AREPO's Config.sh sets FOO, this
script picks it up automatically. Historical species option spellings are
normalized through the compatibility map below.

Usage:
    python3 arepo_config_to_extra_opt.py <path/to/Config.sh> [--src DIR ...]

Prints the compiler flags ready to paste into a `make` invocation, e.g.:

    make -j2 EXTRA_OPT="$(python3 arepo_config_to_extra_opt.py Config.sh)"

Notes / deliberate limitations:
  - Only options that ONEZONE's source actually tests for are emitted.
    Config.sh is full of AREPO-only physics (gravity solver, mesh, sink
    particles, MPI, ...) that has no meaning to the chemistry module and
    must not be blindly passed through -- the whitelist is exactly what
    keeps this safe.
  - A Config.sh option with no value (e.g. "HELIUM") becomes "-DCHEMISTRY_SPECIES_HELIUM".
    An option with a value (e.g. "JEANS_LENGTH_TMAX=1.e4") becomes
    "-DJEANS_LENGTH_TMAX=1.e4", value passed through verbatim.
  - This intentionally says nothing about *which* Config.sh is the right
    one to use for a given comparison -- that is a physics/provenance
    question for a human to answer, not something this script can infer.
"""
import argparse
import re
import sys
from pathlib import Path

DEFAULT_SRC_DIRS = ["src/chemistry", "src/radiation", "src/standalone"]

SPECIES_OPTION_ALIASES = {
    "HELIUM": "CHEMISTRY_SPECIES_HELIUM",
    "HD": "CHEMISTRY_SPECIES_DEUTERIUM",
    "H2_PLUS": "CHEMISTRY_SPECIES_H2_PLUS",
    "HD_PLUS": "CHEMISTRY_SPECIES_HD_PLUS",
    "DMINUS": "CHEMISTRY_SPECIES_D_MINUS",
    "HD_MINUS": "CHEMISTRY_SPECIES_D_MINUS",
    "CARBON_OXYGEN_CHEMISTRY": "CHEMISTRY_SPECIES_CARBON_OXYGEN",
    "CARBON_OXYGEN_HIGH_IONIZATION": "CHEMISTRY_SPECIES_CARBON_OXYGEN_HIGH_IONIZATION",
}

PHYSICS_OPTION_ALIASES = {
    "GAMMA_NEW": "CHEMISTRY_H2_ROTATIONAL_GAMMA",
    "THREEBODY": "CHEMISTRY_THREE_BODY",
    "H2_FORMATION_ON_DUST": "CHEMISTRY_H2_FORMATION_ON_DUST",
    "CASE_B": "CHEMISTRY_CASE_B",
    "MARTIN97": "CHEMISTRY_MARTIN97",
    "TBR_ABEL02": "CHEMISTRY_THREE_BODY_ABEL02",
    "TBR_FLOW07": "CHEMISTRY_THREE_BODY_FLOW07",
    "TBR_FORR13": "CHEMISTRY_THREE_BODY_FORR13",
    "ENZO_CHEM": "CHEMISTRY_ENZO",
    "O05_CHEM": "CHEMISTRY_O05",
    "CHEM_HIRANO": "CHEMISTRY_HIRANO",
    "DUST_COOLING": "CHEMISTRY_DUST_COOLING",
    "MOLECULAR_LINE_COOLING": "CHEMISTRY_MOLECULAR_LINE_COOLING",
    "METAL_LINE_COOLING_SD93": "CHEMISTRY_METAL_LINE_COOLING_SD93",
    "HELIUM_EQ_COOLING": "CHEMISTRY_HELIUM_EQ_COOLING",
    "HELIUM_COOLING": "CHEMISTRY_HELIUM_COOLING",
    "HLINE_COOLING_THIN": "CHEMISTRY_HLINE_COOLING_THIN",
    "CONTINUUM_COOLING_MATSUKOBA": "CHEMISTRY_CONTINUUM_COOLING_MATSUKOBA",
    "CONT_OPACITY": "CHEMISTRY_CONT_OPACITY",
    "H2_OPACITY": "CHEMISTRY_H2_OPACITY",
    "H2_OPACITY_FIT": "CHEMISTRY_H2_OPACITY_FIT",
    "H2_COOLING_G15": "CHEMISTRY_H2_COOLING_G15",
    "H2_COOLING_GP12": "CHEMISTRY_H2_COOLING_GP12",
    "PARTITION_H2_EBM1989": "CHEMISTRY_H2_PARTITION_EBM1989",
    "PARTITION_H2_I87": "CHEMISTRY_H2_PARTITION_I87",
    "PARTITION_H2_R85": "CHEMISTRY_H2_PARTITION_R85",
}

CHEMISTRY_OPTION_ALIASES = {
    "LYMANWERNER": "CHEMISTRY_LW_RADIATION",
    "BACKGROUND_LW_RADIATION": "CHEMISTRY_BACKGROUND_LW_RADIATION",
    "SHIELD_LW": "CHEMISTRY_LW_SELF_SHIELDING",
    "IONIZATION_PHOTON": "IONIZING_PHOTON",

    "CHEM_DIFF_ANALYTIC": "CHEMISTRY_DIFF_ANALYTIC",
    "CHEM_DIFF_NUMERICAL": "CHEMISTRY_DIFF_NUMERICAL",
    "CHEM_DISABLE_CARBON_OXYGEN_CONSERVATION": "CHEMISTRY_DISABLE_CARBON_OXYGEN_CONSERVATION",
    "CHEM_ENERGY_BRACKET_FALLBACK": "CHEMISTRY_ENERGY_BRACKET_FALLBACK",
    "CHEM_EQUILIBRIUM_SPECIES_CONVERGENCE": "CHEMISTRY_EQUILIBRIUM_SPECIES_CONVERGENCE",
    "CHEM_EXPLICIT_SEQUENTIAL": "CHEMISTRY_EXPLICIT_SEQUENTIAL",
    "CHEM_EXPLICIT_SIMULTANEOUS": "CHEMISTRY_EXPLICIT_SIMULTANEOUS",
    "CHEM_FULL_IMPLICIT_EQUILIBRIUM_INITIALIZATION": "CHEMISTRY_FULL_IMPLICIT_EQUILIBRIUM_INITIALIZATION",
    "CHEM_INTERMEDIATE_EQUILIBRIUM": "CHEMISTRY_INTERMEDIATE_EQUILIBRIUM",
    "CHEM_NEWTON_LU_REUSE": "CHEMISTRY_NEWTON_LU_REUSE",
    "CHEM_NEWTON_RESIDUAL_GATE_LIMIT": "CHEMISTRY_NEWTON_RESIDUAL_GATE_LIMIT",
    "CHEM_SAHA_DIAGNOSTICS": "CHEMISTRY_SAHA_DIAGNOSTICS",
    "CHEM_CO_O2_DEPENDENT_DENSITY": "CHEMISTRY_CO_O2_DEPENDENT_DENSITY",
    "ENERGY_IMPLICIT": "CHEMISTRY_ENERGY_IMPLICIT",
    "HYBRID_IMPLICIT": "CHEMISTRY_HYBRID_IMPLICIT",
    "SAHA": "CHEMISTRY_SAHA",
    "INVERSE_REACTIONS": "CHEMISTRY_INVERSE_REACTIONS",
    "CMB_HEATING": "CHEMISTRY_CMB_HEATING",
    "UPDATE_L_CORE_IN_SOLVER": "CHEMISTRY_UPDATE_CORE_LENGTH",
    "JEANS_LENGTH_TMAX": "CHEMISTRY_JEANS_TEMPERATURE_MAX",
    "MAX_ITER_CHEM": "CHEMISTRY_MAX_ITERATIONS",
    "FIX_TEMPERATURE": "CHEMISTRY_FIXED_TEMPERATURE",
    "SET_TMAX": "CHEMISTRY_TEMPERATURE_MAX",
}

OPTION_ALIASES = {**SPECIES_OPTION_ALIASES, **PHYSICS_OPTION_ALIASES, **CHEMISTRY_OPTION_ALIASES}

# Matches: #ifdef FOO / #ifndef FOO / #if defined(FOO) / #elif defined FOO
# (with or without parentheses, arbitrary whitespace).
IFDEF_RE = re.compile(r"^\s*#\s*(?:ifdef|ifndef)\s+(\w+)", re.MULTILINE)
DEFINED_RE = re.compile(r"\bdefined\s*\(?\s*(\w+)\s*\)?")

# Config.sh line: NAME or NAME=VALUE, with an optional trailing "# comment".
# AREPO Config.sh values can contain dots/letters/digits (e.g. "1.e4"),
# so keep the value capture permissive but stop at whitespace or '#'.
CONFIG_LINE_RE = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*(?:=\s*(\S+))?\s*(?:#.*)?$")


def collect_whitelist(repo_root: Path, src_dirs):
    names = set()
    for rel in src_dirs:
        base = repo_root / rel
        if not base.exists():
            continue
        for path in base.rglob("*"):
            if "tests" in path.relative_to(base).parts:
                continue
            if path.suffix not in (".c", ".cc", ".cpp", ".h", ".hpp"):
                continue
            try:
                text = path.read_text(errors="ignore")
            except OSError:
                continue
            text = re.sub(r"\\\r?\n", " ", text)
            text = re.sub(r"/\*.*?\*/|//[^\n]*", "", text, flags=re.DOTALL)
            names.update(IFDEF_RE.findall(text))
            for expression in re.findall(r"^\s*#\s*(?:if|elif)\s+([^\n]+)", text, re.MULTILINE):
                names.update(DEFINED_RE.findall(expression))
                # Also support numeric expressions such as #if OPTION or #elif GA08.
                names.update(re.findall(r"\b[A-Z_][A-Z_0-9]*\b", expression))
    return names


def parse_config_sh(path: Path):
    """Return {NAME: value_or_None} for every active (uncommented) option."""
    active = {}
    for raw_line in path.read_text(errors="ignore").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        m = CONFIG_LINE_RE.match(line)
        if not m:
            continue
        name, value = m.group(1), m.group(2)
        active[name] = value
    return active


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config_sh", type=Path, help="path to an AREPO Config.sh")
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parent.parent,
                        help="ONEZONE repo root (default: parent of this script's directory)")
    parser.add_argument("--src", action="append", dest="src_dirs", default=None,
                        help="source dir (relative to --repo-root) to scan for #ifdef "
                             "whitelist; may be given multiple times "
                             "(default: %s)" % ", ".join(DEFAULT_SRC_DIRS))
    parser.add_argument("--show-skipped", action="store_true",
                        help="also print (to stderr) the Config.sh options that were "
                             "active but not recognized by any ONEZONE #ifdef, for "
                             "auditing")
    args = parser.parse_args()

    src_dirs = args.src_dirs or DEFAULT_SRC_DIRS
    whitelist = collect_whitelist(args.repo_root, src_dirs)
    if not whitelist:
        print("error: found no #ifdef/#ifndef/defined(...) identifiers under %s "
              "(wrong --repo-root?)" % ", ".join(src_dirs), file=sys.stderr)
        sys.exit(2)

    active = parse_config_sh(args.config_sh)

    matched = {}
    for name, value in active.items():
        canonical = OPTION_ALIASES.get(name, name)
        if canonical not in whitelist:
            continue
        if canonical in matched and matched[canonical] != value:
            parser.error("conflicting values for %s" % canonical)
        matched[canonical] = value

    flags = []
    for name in sorted(matched):
        value = matched[name]
        flags.append("-D%s=%s" % (name, value) if value else "-D%s" % name)

    print(" ".join(flags))

    if args.show_skipped:
        skipped = sorted(name for name in active if OPTION_ALIASES.get(name, name) not in matched)
        print("# %d option(s) active in %s but not referenced by ONEZONE source "
              "under %s (skipped):" % (len(skipped), args.config_sh, ", ".join(src_dirs)),
              file=sys.stderr)
        for name in skipped:
            print("#   %s" % name, file=sys.stderr)


if __name__ == "__main__":
    main()
