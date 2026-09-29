"""Configuration regressions: C++ discovery and actual solver selection."""
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("options", ROOT / "tools/arepo_config_to_extra_opt.py")
options = importlib.util.module_from_spec(spec)
sys.dont_write_bytecode = True
spec.loader.exec_module(options)


class BuildOptionsTest(unittest.TestCase):
    def test_cpp_discovery_and_test_exclusion(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / "src/tests").mkdir(parents=True)
            (root / "src/code.cc").write_text("#ifdef CPP_OPTION\n#endif\n#if NUMERIC_OPTION > 0\n#endif\n")
            (root / "src/model.hpp").write_text("#if defined(HEADER_OPTION) && \\\n defined(CONTINUED_OPTION)\n#endif\n/* #ifdef COMMENT_ONLY */\n")
            (root / "src/tests/test.cc").write_text("#ifdef TEST_ONLY\n#endif\n")
            self.assertEqual(options.collect_whitelist(root, ["src"]),
                             {"CPP_OPTION", "NUMERIC_OPTION", "HEADER_OPTION", "CONTINUED_OPTION"})

    def test_actual_cpp_only_option(self):
        whitelist = options.collect_whitelist(ROOT, options.DEFAULT_SRC_DIRS)
        self.assertIn("CHEMISTRY_EXPLICIT_SEQUENTIAL", whitelist)
        self.assertIn("GA08", whitelist)
        self.assertNotIn("EXPECT_HIGH_IONIZATION", whitelist)
        self.assertNotIn("CHEMISTRY_SC", whitelist)

    def test_make_selection(self):
        for mode in [None, "SEQUENTIAL", "SIMULTANEOUS"]:
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as folder:
                cmd = ["make", "-n", "-B", "BUILDDIR=" + folder,
                       folder + "/src/standalone/main.o"]
                if mode:
                    cmd.append("EXPLICIT_CHEMISTRY_UPDATE=" + mode)
                flags = subprocess.check_output(cmd, cwd=ROOT, text=True).split()
                sequential = mode == "SEQUENTIAL"
                self.assertEqual("-DCHEMISTRY_EXPLICIT_SEQUENTIAL" in flags, sequential)
                self.assertEqual("-DCHEMISTRY_EXPLICIT_SIMULTANEOUS" in flags, not sequential)
                self.assertIn("-DCHEMISTRY_ENERGY_BRACKET_FALLBACK", flags)
                for obsolete in ["SINK", "SPLIT", "CMB_HEAT_PROPER", "VARIABLE_GAMMA", "CHEMISTRY_SC"]:
                    self.assertNotIn("-D" + obsolete, flags)

    def test_species_alias_layout(self):
        command = ["c++", "-std=c++17", "-x", "c++", "-E", "-P",
                   "-I" + str(ROOT / "src/chemistry"), "-include", "species.h",
                   "-DMETALLICITY=-4", "/dev/null"]
        old = [name for name in options.SPECIES_OPTION_ALIASES if name != "HD_MINUS"]
        new = sorted(set(options.SPECIES_OPTION_ALIASES.values()))
        old_layout = subprocess.check_output(command + ["-D" + name for name in old], text=True)
        new_layout = subprocess.check_output(command + ["-D" + name for name in new], text=True)
        self.assertEqual(old_layout, new_layout)
        old[old.index("DMINUS")] = "HD_MINUS"
        self.assertEqual(subprocess.check_output(command + ["-D" + name for name in old], text=True), new_layout)
        for flags in [["CHEMISTRY_SPECIES_HD_PLUS"], ["CHEMISTRY_SPECIES_D_MINUS"],
                      ["CHEMISTRY_SPECIES_CARBON_OXYGEN_HIGH_IONIZATION"]]:
            result = subprocess.run(command + ["-D" + flag for flag in flags], capture_output=True)
            self.assertNotEqual(result.returncode, 0)

    def test_species_config_import(self):
        with tempfile.TemporaryDirectory() as folder:
            config = Path(folder) / "Config.sh"
            config.write_text("\n".join(options.SPECIES_OPTION_ALIASES) + "\nMETALLICITY=-4\n")
            flags = subprocess.check_output([sys.executable, str(ROOT / "tools/arepo_config_to_extra_opt.py"),
                                             str(config)], text=True).split()
            expected = {"-D" + name for name in options.SPECIES_OPTION_ALIASES.values()}
            self.assertEqual(set(flags), expected | {"-DMETALLICITY=-4"})
            self.assertEqual(len(flags), len(set(flags)))

    def test_species_make_switch(self):
        for name in ["CHEMISTRY_SPECIES_CARBON_OXYGEN", "CARBON_OXYGEN_CHEMISTRY"]:
            with tempfile.TemporaryDirectory() as folder:
                command = ["make", "-n", "-B", "BUILDDIR=" + folder,
                           folder + "/src/standalone/main.o", name + "=0"]
                flags = subprocess.check_output(command, cwd=ROOT, text=True).split()
                self.assertNotIn("-DCHEMISTRY_SPECIES_CARBON_OXYGEN", flags)
                self.assertNotIn("-DCARBON_OXYGEN_CHEMISTRY", flags)

    def check_physics(self, assertions, flags=(), success=True):
        command = ["c++", "-std=c++17", "-x", "c++", "-fsyntax-only",
                   "-include", str(ROOT / "src/chemistry/chemistry_build_config.h"), "-"]
        result = subprocess.run(command + list(flags), input="\n".join(
            'static_assert(' + expr + ', "configuration");' for expr in assertions),
            text=True, capture_output=True)
        self.assertEqual(result.returncode == 0, success, result.stderr)

    def test_physics_defaults_and_overrides(self):
        self.check_physics([
            "CHEMISTRY_THREE_BODY == 1", "CHEMISTRY_THREE_BODY_FORR13 == 1",
            "CHEMISTRY_H2_COOLING_G15 == 1", "CHEMISTRY_DUST_COOLING == 0",
            "CHEMISTRY_HELIUM_COOLING == 0", "CHEMISTRY_MOLECULAR_LINE_COOLING == 0",
            "CHEMISTRY_H2_PARTITION_EBM1989 == 0", "CHEMISTRY_H2_PARTITION_I87 == 0",
            "CHEMISTRY_H2_PARTITION_R85 == 0"])
        self.check_physics(["CHEMISTRY_HELIUM_COOLING == 1",
                            "CHEMISTRY_MOLECULAR_LINE_COOLING == 1"],
                           ["-DCHEMISTRY_SPECIES_HELIUM", "-DCHEMISTRY_SPECIES_CARBON_OXYGEN", "-DMETALLICITY=-4"])
        self.check_physics(["CHEMISTRY_HELIUM_COOLING == 0"],
                           ["-DCHEMISTRY_SPECIES_HELIUM", "-DCHEMISTRY_HELIUM_COOLING=0"])
        self.check_physics(["CHEMISTRY_H2_FORMATION_ON_DUST == 0", "CHEMISTRY_H2_OPACITY_FIT == 0"],
                           ["-DCHEMISTRY_DUST_COOLING=0", "-DCHEMISTRY_H2_OPACITY=0"])
        self.check_physics(["CHEMISTRY_THREE_BODY_FORR13 == 0", "CHEMISTRY_H2_COOLING_G15 == 0"],
                           ["-DCHEMISTRY_THREE_BODY_ABEL02", "-DCHEMISTRY_H2_COOLING_GP12"])
        self.check_physics(["CHEMISTRY_O05 == 1"], ["-DCHEMISTRY_HIRANO"])
        for flags in [["-DCHEMISTRY_THREE_BODY_ABEL02", "-DCHEMISTRY_THREE_BODY_FLOW07"],
                      ["-DCHEMISTRY_H2_COOLING_GP12", "-DCHEMISTRY_H2_COOLING_G15"],
                      ["-DCHEMISTRY_H2_PARTITION_I87", "-DCHEMISTRY_H2_PARTITION_R85"],
                      ["-DCHEMISTRY_DUST_COOLING=2"],
                      ["-DCHEMISTRY_HELIUM_COOLING=1"],
                      ["-DCHEMISTRY_H2_OPACITY_FIT=0"],
                      ["-DCHEMISTRY_DUST_COOLING=0", "-DCHEMISTRY_H2_FORMATION_ON_DUST=1"]]:
            with self.subTest(flags=flags):
                self.check_physics([], flags, success=False)

    def test_coupled_defaults(self):
        self.check_physics(["CHEMISTRY_DUST_COOLING == 1"], ["-DMETALLICITY=0"])
        self.check_physics(["CHEMISTRY_DUST_COOLING == 0"],
                           ["-DMETALLICITY=-4", "-DCHEMISTRY_DUST_COOLING=0"])
        self.check_physics([], ["-DCHEMISTRY_SPECIES_CARBON_OXYGEN"], success=False)
        self.check_physics(["CHEMISTRY_UPDATE_CORE_LENGTH == 1", "CHEMISTRY_H2_ROTATIONAL_GAMMA == 1"])
        self.check_physics(["CHEMISTRY_UPDATE_CORE_LENGTH == 0", "CHEMISTRY_H2_ROTATIONAL_GAMMA == 0"],
                           ["-DCHEMISTRY_UPDATE_CORE_LENGTH=0", "-DCHEMISTRY_H2_ROTATIONAL_GAMMA=0", "-DGAMMA_NEW"])

    def test_physics_make_overrides(self):
        for variable in ["MOLECULAR_LINE_COOLING", "CHEMISTRY_MOLECULAR_LINE_COOLING"]:
            with tempfile.TemporaryDirectory() as folder:
                flags = subprocess.check_output([
                    "make", "-n", "-B", "BUILDDIR=" + folder,
                    folder + "/src/standalone/main.o", variable + "=0",
                    "EXTRA_OPT=-DCHEMISTRY_THREE_BODY_FLOW07 -DCHEMISTRY_H2_COOLING_GP12"],
                    cwd=ROOT, text=True).split()
                self.assertIn("-DCHEMISTRY_MOLECULAR_LINE_COOLING=0", flags)
                self.assertIn("-DCHEMISTRY_THREE_BODY_FLOW07", flags)
                self.assertIn("-DCHEMISTRY_H2_COOLING_GP12", flags)
                for forbidden in ["-DTBR_FORR13", "-DH2_COOLING_G15", "-DMOLECULAR_LINE_COOLING",
                                  "-DCHEMISTRY_THREE_BODY_FORR13", "-DCHEMISTRY_H2_COOLING_G15"]:
                    self.assertNotIn(forbidden, flags)

    def test_physics_compatibility_and_import(self):
        self.check_physics(["CHEMISTRY_THREE_BODY_FLOW07 == 1", "CHEMISTRY_THREE_BODY_FORR13 == 0"],
                           ["-DTBR_FLOW07"])
        self.check_physics(["CHEMISTRY_DUST_COOLING == 0"],
                           ["-DDUST_COOLING", "-DCHEMISTRY_DUST_COOLING=0"])
        with tempfile.TemporaryDirectory() as folder:
            config = Path(folder) / "Config.sh"
            config.write_text("TBR_ABEL02\nH2_COOLING_GP12\nDUST_COOLING=0\n")
            flags = subprocess.check_output([sys.executable, str(ROOT / "tools/arepo_config_to_extra_opt.py"),
                                             str(config)], text=True).split()
            self.assertEqual(set(flags), {"-DCHEMISTRY_THREE_BODY_ABEL02", "-DCHEMISTRY_H2_COOLING_GP12",
                                          "-DCHEMISTRY_DUST_COOLING=0"})
            self.check_physics(["CHEMISTRY_THREE_BODY_FORR13 == 0", "CHEMISTRY_DUST_COOLING == 0"], flags)

    def test_chemistry_names_and_host_boundary(self):
        numeric = {"UPDATE_L_CORE_IN_SOLVER": "0", "MAX_ITER_CHEM": "37",
                   "JEANS_LENGTH_TMAX": "1234", "SET_TMAX": "5000",
                   "CHEM_NEWTON_RESIDUAL_GATE_LIMIT": "250",
                   "CHEM_CO_O2_DEPENDENT_DENSITY": "1.e8"}
        for old, new in options.CHEMISTRY_OPTION_ALIASES.items():
            value = numeric.get(old, "1")
            with self.subTest(option=old):
                self.check_physics([new + " == " + value], ["-D" + old + "=" + value])
        with tempfile.TemporaryDirectory() as folder:
            config = Path(folder) / "Config.sh"
            config.write_text("ENERGY_IMPLICIT\nCHEM_NEWTON_LU_REUSE\nUPDATE_L_CORE_IN_SOLVER=0\n"
                              "RYTR_NCO\nLYMANWERNER\nMETALLICITY=-4\nMETAL_ADVECTION_SC\n")
            flags = subprocess.check_output([sys.executable, str(ROOT / "tools/arepo_config_to_extra_opt.py"),
                                             str(config)], text=True).split()
            self.assertEqual(set(flags), {"-DCHEMISTRY_ENERGY_IMPLICIT", "-DCHEMISTRY_NEWTON_LU_REUSE",
                                          "-DCHEMISTRY_UPDATE_CORE_LENGTH=0", "-DRYTR_NCO", "-DCHEMISTRY_LW_RADIATION",
                                          "-DMETALLICITY=-4", "-DMETAL_ADVECTION_SC"})
        for variable in ["UPDATE_L_CORE", "CHEMISTRY_UPDATE_CORE_LENGTH"]:
            with tempfile.TemporaryDirectory() as folder:
                flags = subprocess.check_output(["make", "-n", "-B", "BUILDDIR=" + folder,
                       folder + "/src/standalone/main.o", variable + "=0"], cwd=ROOT, text=True).split()
                self.assertIn("-DCHEMISTRY_UPDATE_CORE_LENGTH=0", flags)
                self.assertNotIn("-DUPDATE_L_CORE_IN_SOLVER", flags)

    def test_radiation_names_and_local_host_option(self):
        old = ["LYMANWERNER", "BACKGROUND_LW_RADIATION", "SHIELD_LW", "IONIZATION_PHOTON"]
        new = ["CHEMISTRY_LW_RADIATION", "CHEMISTRY_BACKGROUND_LW_RADIATION",
               "CHEMISTRY_LW_SELF_SHIELDING", "IONIZING_PHOTON"]
        command = ["c++", "-std=c++17", "-x", "c++", "-E", "-P",
                   "-I" + str(ROOT / "src/chemistry"), "-include", "chemistry_types.h", "/dev/null"]
        self.assertEqual(subprocess.check_output(command + ["-D" + n for n in old], text=True),
                         subprocess.check_output(command + ["-D" + n for n in new], text=True))
        with tempfile.TemporaryDirectory() as folder:
            config = Path(folder) / "Config.sh"
            config.write_text("\n".join(old + ["LW_RADIATION_LOCAL", "RYTR_NCO"]) + "\n")
            flags = subprocess.check_output([sys.executable, str(ROOT / "tools/arepo_config_to_extra_opt.py"),
                                             str(config)], text=True).split()
            self.assertEqual(set(flags), {"-D" + n for n in new + ["LW_RADIATION_LOCAL", "RYTR_NCO"]})

    def test_aliases_and_conflicts(self):
        header = ROOT / "src/chemistry/chemistry_build_config.h"
        command = ["c++", "-std=c++17", "-x", "c++", "-E", "-dM", "-include", str(header), "/dev/null"]
        definitions = subprocess.check_output(command + ["-DFLAG_BRENT", "-DDEBUG_CHEMISTRY_SC"], text=True)
        self.assertIn("#define CHEMISTRY_ENERGY_BRACKET_FALLBACK", definitions)
        self.assertIn("#define CHEMISTRY_SAHA_DIAGNOSTICS", definitions)
        result = subprocess.run(command + ["-DCHEMISTRY_EXPLICIT_SEQUENTIAL", "-DCHEMISTRY_EXPLICIT_SIMULTANEOUS"],
                                text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Choose only one explicit chemistry update method", result.stderr)


if __name__ == "__main__":
    unittest.main()
