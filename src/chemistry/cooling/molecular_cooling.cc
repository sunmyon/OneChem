#include "chemistry_model.hpp"
#include <new>
#include <algorithm>
#include <initializer_list>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "chemistry_core.h"
#include "chemistry_constants.h"
#include "chemistry_tables.h"
#include "molecular_cooling.h"

#if CHEMISTRY_MOLECULAR_LINE_COOLING

#define MOLECULAR_COOLING_MAX_T 11
#define MOLECULAR_COOLING_MAX_COLUMN 11

class MolecularCoolingTable final {
public:
  MolecularCoolingTable() = default;
  MolecularCoolingTable(const char *path, int temperatures, int columns,
                        double column_min, double column_step,
                        std::initializer_list<double> nodes,
                        std::initializer_list<double> low_density)
      : filename(path), temperature_count(temperatures), column_count(columns),
        log_column_min(column_min), log_column_step(column_step) {
    std::copy(nodes.begin(), nodes.end(), temperature);
    std::copy(low_density.begin(), low_density.end(), log_L0);
  }
  int read_molecular_cooling_table(const char *table_root);
  double net_cooling_per_volume(double molecule_density, double collider_density, double temperature, double radiation_temperature, double effective_column_value) const;
private:
  int read_table_block(FILE *file, double values[MOLECULAR_COOLING_MAX_COLUMN][MOLECULAR_COOLING_MAX_T], int change_sign);
  void interpolation_coordinate(double temperature, double effective_column, int *temperature_index, int *column_index, double *temperature_fraction, double *column_fraction) const;
  static double interpolate_temperature(const double *values, int temperature_index, double temperature_fraction);
  static double interpolate_table(const double values[MOLECULAR_COOLING_MAX_COLUMN][MOLECULAR_COOLING_MAX_T], int temperature_index, int column_index, double temperature_fraction, double column_fraction);
  double cooling_per_molecule(double temperature, double effective_column, double collider_density) const;
  const char *filename{};
  int temperature_count{};
  int column_count{};
  double log_column_min{};
  double log_column_step{};
  double temperature[MOLECULAR_COOLING_MAX_T]{};
  double log_L0[MOLECULAR_COOLING_MAX_T]{};
  double log_L_LTE[MOLECULAR_COOLING_MAX_COLUMN][MOLECULAR_COOLING_MAX_T]{};
  double log_nhalf[MOLECULAR_COOLING_MAX_COLUMN][MOLECULAR_COOLING_MAX_T]{};
  double alpha[MOLECULAR_COOLING_MAX_COLUMN][MOLECULAR_COOLING_MAX_T]{};
};

struct ChemistryCoolingTables {
public:
  ChemistryCoolingTables() = default;
  ChemistryCoolingTables(const ChemistryCoolingTables &) = delete;
  ChemistryCoolingTables &operator=(const ChemistryCoolingTables &) = delete;
  ~ChemistryCoolingTables();

  static chemistry::ResourceOwner<ChemistryCoolingTables> create(const char *table_root) noexcept;


  void evaluate_molecular_line_cooling(double temperature,
                                     double radiation_temperature,
                                     double length, double n_nuclei,
                                     const double *abundance,
                                     double *cooling_CO, double *cooling_OH,
                                     double *cooling_H2O) const;

private:

  static double effective_column(double number_density, double length, double temperature, double molecular_weight);

  MolecularCoolingTable CO;
  MolecularCoolingTable OH;
  MolecularCoolingTable H2O;

};

int MolecularCoolingTable::read_table_block(FILE *file, double values[MOLECULAR_COOLING_MAX_COLUMN][MOLECULAR_COOLING_MAX_T], int change_sign) {
  MolecularCoolingTable *table = this;
  for(int column = 0; column < table->column_count; column++)
    for(int temperature = 0; temperature < table->temperature_count; temperature++) {
      if(fscanf(file, "%lf", &values[column][temperature]) != 1)
        return 0;
      if(change_sign)
        values[column][temperature] = -values[column][temperature];
    }
  return 1;
}

int MolecularCoolingTable::read_molecular_cooling_table(const char *table_root) {
  MolecularCoolingTable *table = this;
  FILE *file = open_chemistry_table(table_root, table->filename);
  if(file == NULL)
    return 0;

  int valid = read_table_block(file, table->log_L_LTE, 1) &&
              read_table_block(file, table->log_nhalf, 0) &&
              read_table_block(file, table->alpha, 0);
  double extra;
  if(valid && fscanf(file, "%lf", &extra) == 1)
    valid = 0;
  fclose(file);

  if(!valid)
    fprintf(stderr, "Invalid molecular cooling table: %s\n", table->filename);
  return valid;
}

chemistry::ResourceOwner<ChemistryCoolingTables> ChemistryCoolingTables::create(const char *table_root) noexcept {

  chemistry::ResourceOwner<ChemistryCoolingTables> tables{new(std::nothrow) ChemistryCoolingTables{}};
  if(tables == NULL)
    return {};

  tables->CO = MolecularCoolingTable {
    "COline_cooling.dat",
    11,
    11,
    14.0,
    0.5,
    {10., 20., 30., 50., 80., 100., 300., 600., 1000., 1500., 2000.},
    {-24.77, -24.38, -24.21, -24.03, -23.89, -23.82, -23.42, -23.13, -22.91, -22.63, -22.28}
  };
  tables->OH = MolecularCoolingTable {
    "OHline_cooling.dat",
    6,
    9,
    10.0,
    1.0,
    {30., 50., 80., 100., 300., 600.},
    {-25.31, -24.53, -24.02, -23.82, -23.16, -22.90}
  };
  tables->H2O = MolecularCoolingTable {
    "H2Oline_cooling.dat",
    6,
    10,
    10.0,
    1.0,
    {100., 200., 400., 1000., 2000., 4000.},
    {-24.35, -23.87, -23.42, -22.88, -22.50, -22.14}
  };

  if(!tables->CO.read_molecular_cooling_table(table_root) ||
     !tables->OH.read_molecular_cooling_table(table_root) ||
     !tables->H2O.read_molecular_cooling_table(table_root)) {

    return {};
  }
  return tables;
}

chemistry::ResourceOwner<ChemistryCoolingTables> chemistry::load_molecular_cooling(const char *root) noexcept {
  return ChemistryCoolingTables::create(root);
}


ChemistryCoolingTables::~ChemistryCoolingTables() {

}

void chemistry::ModelResourceDeleter::operator()(ChemistryCoolingTables *p) const noexcept { delete p; }




void MolecularCoolingTable::interpolation_coordinate(double temperature, double effective_column, int *temperature_index, int *column_index, double *temperature_fraction, double *column_fraction) const {
  const MolecularCoolingTable *table = this;
  double log_temperature = log10(fmax(temperature, table->temperature[0]));
  double log_temperature_max = log10(table->temperature[table->temperature_count - 1]);
  if(log_temperature > log_temperature_max)
    log_temperature = log_temperature_max;

  int lower = 0;
  int upper = table->temperature_count - 1;
  while(upper - lower > 1) {
    int middle = (lower + upper) / 2;
    if(log_temperature < log10(table->temperature[middle]))
      upper = middle;
    else
      lower = middle;
  }
  *temperature_index = lower;
  double log_temperature_lower = log10(table->temperature[lower]);
  double log_temperature_upper = log10(table->temperature[lower + 1]);
  *temperature_fraction = (log_temperature - log_temperature_lower) /
                          (log_temperature_upper - log_temperature_lower);

  double log_column = log10(fmax(effective_column, pow(10.0, table->log_column_min)));
  double log_column_max = table->log_column_min + table->log_column_step * (table->column_count - 1);
  if(log_column > log_column_max)
    log_column = log_column_max;
  double column_position = (log_column - table->log_column_min) / table->log_column_step;
  int column = (int)floor(column_position);
  if(column >= table->column_count - 1)
    column = table->column_count - 2;
  *column_index = column;
  *column_fraction = column_position - column;
}

double MolecularCoolingTable::interpolate_temperature(const double *values, int temperature_index, double temperature_fraction) {
  return values[temperature_index] + temperature_fraction *
         (values[temperature_index + 1] - values[temperature_index]);
}

double MolecularCoolingTable::interpolate_table(const double values[MOLECULAR_COOLING_MAX_COLUMN][MOLECULAR_COOLING_MAX_T], int temperature_index, int column_index, double temperature_fraction, double column_fraction) {
  double lower = interpolate_temperature(values[column_index], temperature_index, temperature_fraction);
  double upper = interpolate_temperature(values[column_index + 1], temperature_index, temperature_fraction);
  return lower + column_fraction * (upper - lower);
}

double MolecularCoolingTable::cooling_per_molecule(double temperature, double effective_column, double collider_density) const {
  const MolecularCoolingTable *table = this;
  if(temperature <= 0.0 || collider_density <= 0.0)
    return 0.0;

  int temperature_index, column_index;
  double temperature_fraction, column_fraction;
  interpolation_coordinate(temperature, effective_column, &temperature_index, &column_index,
                           &temperature_fraction, &column_fraction);

  double log_L0 = interpolate_temperature(table->log_L0, temperature_index, temperature_fraction);
  double log_L_LTE = interpolate_table(table->log_L_LTE, temperature_index, column_index,
                                       temperature_fraction, column_fraction);
  double log_nhalf = interpolate_table(table->log_nhalf, temperature_index, column_index,
                                       temperature_fraction, column_fraction);
  double alpha = interpolate_table(table->alpha, temperature_index, column_index,
                                   temperature_fraction, column_fraction);

  double L0 = pow(10.0, log_L0);
  double L_LTE = pow(10.0, log_L_LTE);
  double nhalf = pow(10.0, log_nhalf);
  double density_ratio = collider_density / nhalf;
  double inverse_coefficient = 1.0 / L0 + collider_density / L_LTE +
    pow(density_ratio, alpha) / L0 * (1.0 - nhalf * L0 / L_LTE);

  if(!isfinite(inverse_coefficient) || inverse_coefficient <= 0.0)
    return 0.0;
  return collider_density / inverse_coefficient;
}

double ChemistryCoolingTables::effective_column(double number_density, double length, double temperature, double molecular_weight) {
  if(number_density <= 0.0 || length <= 0.0 || temperature <= 0.0)
    return 0.0;
  double doppler_width_km_s = sqrt(2.0 * BOLTZMANN * temperature /
                                  (molecular_weight * PROTONMASS)) / 1.e5;
  return number_density * length / doppler_width_km_s;
}

double MolecularCoolingTable::net_cooling_per_volume(double molecule_density, double collider_density, double temperature, double radiation_temperature, double effective_column_value) const {
  const MolecularCoolingTable *table = this;
  if(molecule_density <= 0.0)
    return 0.0;
  double cooling = cooling_per_molecule(temperature, effective_column_value, collider_density);
  if(radiation_temperature > 0.0)
    cooling -= cooling_per_molecule(radiation_temperature, effective_column_value, collider_density);
  return molecule_density * cooling;
}

void ChemistryCoolingTables::evaluate_molecular_line_cooling(double temperature,
                                     double radiation_temperature,
                                     double length, double n_nuclei,
                                     const double *abundance,
                                     double *cooling_CO, double *cooling_OH,
                                     double *cooling_H2O) const {
  const ChemistryCoolingTables *tables = this;
#ifdef CHEMISTRY_PROFILE
  double profile_timer = profile_start();
#endif
  double nH2 = n_nuclei * abundance[sH2I];
  double nCO = n_nuclei * abundance[sCO];
  double nOH = n_nuclei * abundance[sOH];
  double nH2O = n_nuclei * abundance[sH2O];

  double column_CO = effective_column(nCO, length, temperature, 28.0);
  double column_OH = effective_column(nOH, length, temperature, 17.0);
  double column_H2O = effective_column(nH2O, length, temperature, 18.0);

  *cooling_CO = tables->CO.net_cooling_per_volume(nCO, nH2,
      temperature, radiation_temperature, column_CO);
  *cooling_OH = tables->OH.net_cooling_per_volume(nOH, nH2,
      temperature, radiation_temperature, column_OH);
  *cooling_H2O = tables->H2O.net_cooling_per_volume(nH2O, nH2,
      temperature, radiation_temperature, column_H2O);
#ifdef CHEMISTRY_PROFILE
  profile_end(CHEM_PROFILE_MOLECULAR_COOLING, profile_timer);
#endif
}

void evaluate_molecular_line_cooling(const ChemistryCoolingTables *tables,
                                     double temperature,
                                     double radiation_temperature,
                                     double length, double n_nuclei,
                                     const double *abundance,
                                     double *cooling_CO, double *cooling_OH,
                                     double *cooling_H2O) { tables->evaluate_molecular_line_cooling(temperature, radiation_temperature, length, n_nuclei, abundance, cooling_CO, cooling_OH, cooling_H2O); }

#endif
