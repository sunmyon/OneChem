#!/bin/sh
# Compile C and C++ test inputs separately, then link with the C++ runtime.
set -eu
build=$(mktemp -d "${TMPDIR:-/tmp}/chemistry-core-test.XXXXXX")
trap 'rm -rf "$build"' EXIT HUP INT TERM
index=0
for source in "$@"; do
  case "$source" in
    *.cc) compiler=${CXX:-c++}; standard=-std=c++17 ;;
    *.c) compiler=${CC:-cc}; standard=-std=c11 ;;
    *) echo "Unsupported test source: $source" >&2; exit 1 ;;
  esac
  $compiler $standard ${CPPFLAGS:-} ${CHEM_TEST_FLAGS:-} \
    -c "$source" -o "$build/$index.o"
  index=$((index + 1))
done
${CXX:-c++} "$build"/*.o ${CHEM_TEST_LIBS:-} -o "$build/test"
"$build/test"
