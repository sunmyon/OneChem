#include "chemistry_model.hpp"
#include "include/chemistry/chemistry.h"
#include "tables/chemistry_tables.h"
#include "cooling/cooling_func.h"
#include "chemistry_constants.h"
#include "cooling/dust_cooling.h"
#include "reaction_network/reaction_network.h"

#include <cstdio>
#include <cmath>
#include <limits>
#include <new>
#include <string>
#include <type_traits>
#include <utility>

// Inject failures into each nothrow Model/resource/opaque-handle allocation.
static int allocations_before_failure = -1;
void *operator new(std::size_t size, const std::nothrow_t &) noexcept
{
  if(allocations_before_failure == 0)
    return nullptr;
  if(allocations_before_failure > 0)
    --allocations_before_failure;
  try { return ::operator new(size); }
  catch(...) { return nullptr; }
}

void *operator new[](std::size_t size, const std::nothrow_t &) noexcept
{
  if(allocations_before_failure == 0) return nullptr;
  if(allocations_before_failure > 0) --allocations_before_failure;
  try { return ::operator new[](size); }
  catch(...) { return nullptr; }
}

using chemistry::Model;
static_assert(!std::is_copy_constructible<Model>::value, "unique table owner");
static_assert(!std::is_move_constructible<Model>::value, "stable borrowed address");
static_assert(std::is_nothrow_move_constructible<Model::Handle>::value, "movable handle");
static_assert(std::is_const<std::remove_pointer_t<decltype(ChemRunView::rate_tables)>>::value,
              "evaluation cannot modify rate tables");
static_assert(std::is_const<std::remove_pointer_t<decltype(ChemRunView::network)>>::value,
              "evaluation cannot modify the network");

#define CHECK(condition) do { \
  if(!(condition)) { \
    std::fprintf(stderr, "Model lifecycle failed: %s\n", #condition); \
    return 1; \
  } \
} while(0)

int main(int argc, char **argv)
{
  CHECK(argc == 4);
#if CHEMISTRY_CONTINUUM_COOLING_MATSUKOBA && CHEMISTRY_HELIUM_COOLING
  // Isolate the He recombination term and the H-minus terms independently:
  // enabling helium must not overwrite the H-minus fit coefficients.
  for(double temperature : {1000., 10000., 100000.}) {
    double abund[N_SP]{};
    const double logT = std::log10(temperature);
    abund[sELEC] = 0.1;
    abund[sHeII] = 0.05;
    double expected = 4 * M_PI * std::pow(10., -21.4171 - 0.995727 * logT +
                                         0.0714639 * logT * logT) * 0.1 * 0.05;
    double actual = continuum_cooling_rate_thin(temperature, PROTONMASS, abund, 0);
    CHECK(std::fabs(actual / expected - 1.) < 1.e-12);
    abund[sHeII] = 0.;
    abund[sHI] = 0.8;
    expected = 4 * M_PI * 0.1 * 0.8 *
      (std::pow(10., -31.211 + 1.1678 * logT - 0.030473 * logT * logT) +
       std::pow(10., -38.564 + 3.4681 * logT - 0.17880 * logT * logT));
    actual = continuum_cooling_rate_thin(temperature, PROTONMASS, abund, 0);
    CHECK(std::fabs(actual / expected - 1.) < 1.e-12);
  }
#endif

  // Linear values on log-spaced nodes: interpolate in T, not log(T).
  const double nodes[] = {1., 10., 100.};
  const double values[] = {3., 21., 201.};
  const chemistry::LogGridView grid(nodes, 3, std::log(1.), std::log(100.),
                                     std::log(10.));
  for(double temperature : {0.1, 1., 5., 10., 55., 100., 200.}) {
    const auto sample = grid.locate(temperature);
    CHECK(sample.lower >= 0 && sample.lower < 2);
    CHECK(std::fabs(sample.interpolate(values) - (2. * temperature + 1.)) < 1.e-12);
  }
  CHECK(grid.locate(0.1).lower == 0 && grid.locate(0.1).fraction < 0.);
  CHECK(grid.locate(200.).lower == 1 && grid.locate(200.).fraction > 1.);
  for(double temperature : {std::nextafter(10., 0.), std::nextafter(10., 100.)}) {
    // libm rounding can choose either adjacent interval at a node.
    CHECK(std::fabs(grid.locate(temperature).interpolate(values) -
                    (2. * temperature + 1.)) < 1.e-12);
  }

  ChemistryConfig config{};
  chemistry_config_defaults(&config);
  std::string root = argv[1];
  config.table_root = root.c_str();
  CHECK(config.print_configuration == 1);
  auto first = Model::create(config);
  CHECK(first.status == CHEMISTRY_STATUS_OK && first.model);
  config.print_configuration = 0; // Subsequent successful creations must be silent.
  const auto *view = &first.model->run();
  CHECK(first.model->network() == view->network);
  CHECK(first.model->rates().data() == view->rate_tables);
  CHECK(first.model->molecular_cooling() == view->cooling_tables);
  CHECK(first.model->dust_opacity() == view->dust_tables);
  CHECK(first.model->continuum_opacity() == view->continuum_opacity_tables);
  // The alternate root also uses historical species names in the reaction table.
  auto alternate_config = config;
  alternate_config.table_root = argv[3];
  auto alternate = Model::create(alternate_config);
  CHECK(alternate.status == CHEMISTRY_STATUS_OK && alternate.model);
  CHECK(alternate.model->network()->size() == first.model->network()->size());
#if CHEMISTRY_METAL_LINE_COOLING_SD93
  // Two roots must not overwrite each other's metal cooling data.
  const double metal_rate = get_cooling_rate_SD93(view->continuum_opacity_tables, 15000.);
  CHECK(std::isfinite(metal_rate) && metal_rate > 0.);
  const double alternate_rate = get_cooling_rate_SD93(
      alternate.model->run().continuum_opacity_tables, 15000.);
  CHECK(std::fabs(alternate_rate / metal_rate - 2.) < 1.e-13);
  CHECK(get_cooling_rate_SD93(view->continuum_opacity_tables, 15000.) == metal_rate);
#endif
  alternate.model.reset();
  auto moved = std::move(first.model);
  CHECK(!first.model && &moved->run() == view);

  config.timescale_safety_factor = 0.25;
  auto second = Model::create(config);
  CHECK(second.status == CHEMISTRY_STATUS_OK && second.model);
  CHECK(view->network != second.model->run().network);
  CHECK(view->rate_tables != second.model->run().rate_tables);
  CHECK(view->timescale_safety_factor == 0.1);
  CHECK(second.model->run().timescale_safety_factor == 0.25);
  const double first_temperature = second.model->run().rate_tables->temperature[0];
  const auto reaction_count = (second.model->network())->size();
#if CHEMISTRY_DUST_COOLING
  double rosseland_before, planck_before;
  chemistry_dust_opacities_from_tables(second.model->dust_opacity(), 300.,
                                       &rosseland_before, &planck_before);
#endif
  root.assign(root.size(), '?'); // table_root is not retained after loading.
  moved.reset();
  CHECK(second.model->run().rate_tables->temperature[0] == first_temperature);
  CHECK((second.model->network())->size() == reaction_count);
#if CHEMISTRY_DUST_COOLING
  double rosseland_after, planck_after;
  chemistry_dust_opacities_from_tables(second.model->dust_opacity(), 300.,
                                       &rosseland_after, &planck_after);
  CHECK(rosseland_before == rosseland_after && planck_before == planck_after);
#endif


  config.table_root = argv[1];
  config.timescale_safety_factor = std::numeric_limits<double>::quiet_NaN();
  auto invalid = Model::create(config);
  CHECK(!invalid.model && invalid.status == CHEMISTRY_STATUS_INVALID_ARGUMENT);
  config.timescale_safety_factor = 0.1;

  allocations_before_failure = 0;
  auto no_owner = Model::create(config);
  allocations_before_failure = -1;
  CHECK(!no_owner.model && no_owner.status == CHEMISTRY_STATUS_OUT_OF_MEMORY);

  // Exercise every Model/resource allocation without an opaque C context.
  bool native_success = false;
  for(int fail_after = 1; fail_after < 2 * N_SP + 32; ++fail_after) {
    allocations_before_failure = fail_after;
    auto candidate = Model::create(config);
    allocations_before_failure = -1;
    if(candidate.status == CHEMISTRY_STATUS_OK) {
      CHECK(candidate.model);
      native_success = true;
      break;
    }
    CHECK(!candidate.model);
    CHECK(candidate.status == CHEMISTRY_STATUS_OUT_OF_MEMORY || candidate.status == CHEMISTRY_STATUS_IO_ERROR);
  }
  CHECK(native_success);
  // Missing a late-loaded dust table exercises cleanup after network, rates,
  // molecular cooling and continuum tables have already been constructed.
#if CHEMISTRY_DUST_COOLING
  config.table_root = argv[2];
  for(int attempt = 0; attempt < 3; ++attempt) {
    auto failed = Model::create(config);
    CHECK(!failed.model && failed.status == CHEMISTRY_STATUS_IO_ERROR);
  }
#endif
  config.table_root = argv[1];
  std::puts("PASS: Model ownership, stable views, independent runs and failure recovery");
  return 0;
}
