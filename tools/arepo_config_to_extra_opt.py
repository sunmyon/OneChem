#!/usr/bin/env python3
"""Extract the ONEZONE-relevant compile options from an AREPO Config.sh.

Problem this solves: ONEZONE and AREPO's chemistry_ONECHEM module are meant
to share chemistry_solver.c byte-for-byte, and their behavior is controlled
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
script picks it up automatically, with no hardcoded flag list to keep in
sync by hand.

Usage:
    python3 arepo_config_to_extra_opt.py <path/to/Config.sh> [--src DIR ...]

Prints an EXTRA_OPT='...' line ready to paste into a `make` invocation, e.g.:

    make -j2 EXTRA_OPT="$(python3 arepo_config_to_extra_opt.py Config.sh)"

Notes / deliberate limitations:
  - Only options that ONEZONE's source actually tests for are emitted.
    Config.sh is full of AREPO-only physics (gravity solver, mesh, sink
    particles, MPI, ...) that has no meaning to the chemistry module and
    must not be blindly passed through -- the whitelist is exactly what
    keeps this safe.
  - A Config.sh option with no value (e.g. "HELIUM") becomes "-DHELIUM".
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
            if path.suffix not in (".c", ".h"):
                continue
            try:
                text = path.read_text(errors="ignore")
            except OSError:
                continue
            names.update(IFDEF_RE.findall(text))
            names.update(DEFINED_RE.findall(text))
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

    matched = {name: value for name, value in active.items() if name in whitelist}

    flags = []
    for name in sorted(matched):
        value = matched[name]
        flags.append("-D%s=%s" % (name, value) if value else "-D%s" % name)

    print(" ".join(flags))

    if args.show_skipped:
        skipped = sorted(set(active) - set(matched))
        print("# %d option(s) active in %s but not referenced by ONEZONE source "
              "under %s (skipped):" % (len(skipped), args.config_sh, ", ".join(src_dirs)),
              file=sys.stderr)
        for name in skipped:
            print("#   %s" % name, file=sys.stderr)


if __name__ == "__main__":
    main()
