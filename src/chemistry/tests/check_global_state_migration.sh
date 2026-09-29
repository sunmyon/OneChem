#!/bin/sh
set -eu

network_source=reaction_network/reaction_network.cc

if grep -n -E '^static (ChemReaction \*reactions|size_t reaction_count|ChemSpeciesReactionIndex species_reaction_index)' "$network_source"; then
  echo "error: reaction-network mutable file-scope state was reintroduced" >&2
  exit 1
fi

if grep -R -n -E '(^|[^A-Za-z0-9_])N_RE([^A-Za-z0-9_]|$)' . \
    --include='*.c' --include='*.cc' --include='*.cpp' --include='*.h' --include='*.hpp'; then
  echo "error: legacy global reaction count was reintroduced" >&2
  exit 1
fi

legacy_rate_globals='T_table|k_table|k_table_inverse|ciHIa|ciHeIa|ciHeIIa|ciHeISa|reHIIa|brema|H2_CIE|ceHIa|ceHeIa|ceHeIIa|reHeII1a|reHeII2a|reHeIIIa|chem_e'
if grep -R -n -E "^(extern[[:space:]]+)?double[[:space:]]+(${legacy_rate_globals})(\[|[[:space:]]*[;,=])" . \
    --include='*.c' --include='*.cc' --include='*.cpp' --include='*.h' --include='*.hpp'; then
  echo "error: a legacy mutable rate-table global was reintroduced" >&2
  exit 1
fi

if grep -n -E '^static MolecularCoolingTable (CO_table|OH_table|H2O_table)' \
    cooling/molecular_cooling.cc; then
  echo "error: a molecular cooling table global was reintroduced" >&2
  exit 1
fi

legacy_dust_globals='n_Tdust|Kappa_Tdust|Kappa_Dust_Rossland|Kappa_Dust_Planck|dust_log_T_min|dust_inv_dlog_T'
if grep -n -E "^(static[[:space:]]+)?(int|double|double[[:space:]]*\*)[[:space:]]+(${legacy_dust_globals})(\[|[[:space:]]*[;,=])" \
    cooling/dust_cooling.cc; then
  echo "error: a legacy mutable dust-opacity global was reintroduced" >&2
  exit 1
fi

legacy_continuum_opacity_globals='Kappa_Tgas|Kappa_rho|Kappa_Rossland|Kappa_Planck|RhoMin_tau|RhoMax_tau|TMin_tau|TMax_tau|n_T|n_rho'
if grep -n -E "^(int|double|double[[:space:]]*\*)[[:space:]]+(${legacy_continuum_opacity_globals})(\[|[[:space:]]*[;,=])" \
    cooling/cooling_func.cc; then
  echo "error: a legacy mutable continuum-opacity global was reintroduced" >&2
  exit 1
fi

if grep -R -n -E 'ChemistryRadiationTables|radiation_tables|react_rate_FUV|get_Hminus_PD_rate' . \
    --include='*.c' --include='*.cc' --include='*.cpp' --include='*.h' --include='*.hpp'; then
  echo "error: host background-radiation construction leaked into chemistry" >&2
  exit 1
fi

if grep -R -n -E '^(extern[[:space:]]+)?const double (XH|XHe|XHD)([[:space:]]|=)' . \
    --include='*.c' --include='*.cc' --include='*.cpp' --include='*.h' --include='*.hpp'; then
  echo "error: chemistry composition was reintroduced as external state" >&2
  exit 1
fi

if grep -R -n -E '(^|[^A-Za-z0-9_])(FlagPrint|count_print)([^A-Za-z0-9_]|$)' . \
    --include='*.c' --include='*.cc' --include='*.cpp' --include='*.h' --include='*.hpp'; then
  echo "error: process-wide diagnostic state was reintroduced" >&2
  exit 1
fi

if grep -R -n -E 'chemistry_table_root|table_root_is_set' . \
    --include='*.c' --include='*.cc' --include='*.cpp' --include='*.h' --include='*.hpp'; then
  echo "error: process-wide table-root state was reintroduced" >&2
  exit 1
fi

if grep -R -n -E 'open_chemistry_table\([^,()]+\)' . \
    --include='*.c' --include='*.cc' --include='*.cpp' --include='*.h' --include='*.hpp'; then
  echo "error: chemistry table opener must receive table_root explicitly" >&2
  exit 1
fi

echo "run tables and diagnostics have explicit owners; composition is immutable ABI metadata"
