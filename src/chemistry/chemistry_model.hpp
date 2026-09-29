#pragma once

#include <memory>
#include "chemistry_core.h"
#include "tables/chemistry_tables.h"
#include "include/chemistry/chemistry_config.h"

namespace chemistry {

// A borrowed, logarithmically spaced axis. count >= 2, positive ordered nodes,
// and a positive log step are preconditions supplied by the table builder.
// Lookup clamps the interval, NOT the input: evaluation extrapolates at either end.
class LogGridView final {
public:
  struct Sample {
    int lower;
    double fraction;
    double interpolate(const double *values) const noexcept {
      return values[lower] + fraction * (values[lower + 1] - values[lower]);
    }
  };
  LogGridView(const double *nodes, int count, double log_min,
              double log_max, double log_step) noexcept;
  Sample locate(double value) const noexcept; // finite, positive input
private:
  const double *nodes_;
  int count_;
  double log_min_, log_max_, log_step_;
};

// Non-owning evaluation view; the Model must outlive it.
// Locate once per temperature and reuse the interval for all reaction channels.
class RateTableView final {
public:
  RateTableView(const ChemistryRateTables &tables, double temperature) noexcept;
  double forward(int reaction) const noexcept {
    return sample_.interpolate(tables_.forward[reaction]);
  }
  double inverse(int reaction) const noexcept {
    return sample_.interpolate(tables_.inverse[reaction]);
  }
  int lower_index() const noexcept { return sample_.lower; }
  double fraction() const noexcept { return sample_.fraction; }
private:
  const ChemistryRateTables &tables_;
  LogGridView::Sample sample_;
};

// Owns ChemistryRateTables; evaluation code receives only a const borrowed pointer.
class RateTables final {
public:
  RateTables() noexcept;
  ~RateTables() noexcept;
  RateTables(const RateTables &) = delete;
  RateTables &operator=(const RateTables &) = delete;
  const ChemistryRateTables *data() const noexcept { return data_.get(); }

private:
  friend class Model;
  std::unique_ptr<ChemistryRateTables> data_;
  ChemistryStatus build() noexcept;
  static int tabulate(ChemistryRateTables *tables);
};

// Deletion is defined beside each complete resource type. Model owns those
// objects directly; no separate per-resource ownership wrapper is needed.
struct ModelResourceDeleter {
  void operator()(ChemReactionNetwork *) const noexcept;
  void operator()(ChemistryCoolingTables *) const noexcept;
  void operator()(ChemistryDustTables *) const noexcept;
  void operator()(ChemistryContinuumOpacityTables *) const noexcept;
};

template<class T>
using ResourceOwner = std::unique_ptr<T, ModelResourceDeleter>;

ResourceOwner<ChemReactionNetwork> load_network(const char *root) noexcept;
#if CHEMISTRY_MOLECULAR_LINE_COOLING
ResourceOwner<ChemistryCoolingTables> load_molecular_cooling(const char *root) noexcept;
#endif
#if CHEMISTRY_DUST_COOLING
ResourceOwner<ChemistryDustTables> load_dust_opacity(const char *root) noexcept;
#endif
#if CHEMISTRY_CONTINUUM_COOLING_MATSUKOBA || CHEMISTRY_METAL_LINE_COOLING_SD93
ResourceOwner<ChemistryContinuumOpacityTables> load_continuum_opacity(const char *root) noexcept;
#endif

// Run-lifetime aggregate: explicit owners, then one borrowed evaluation view.
class Model final {
public:
  using Handle = std::unique_ptr<const Model>;
  struct Creation;
  static Creation create(const ChemistryConfig &config) noexcept;
  ~Model() noexcept;
  Model(const Model &) = delete;
  Model &operator=(const Model &) = delete;
  Model(Model &&) = delete;
  Model &operator=(Model &&) = delete;

  const ChemRunView &run() const noexcept { return view_; }
  const ChemReactionNetwork *network() const noexcept { return network_.get(); }
  const RateTables &rates() const noexcept { return rates_; }
  const ChemistryCoolingTables *molecular_cooling() const noexcept { return molecular_cooling_.get(); }
  const ChemistryDustTables *dust_opacity() const noexcept { return dust_opacity_.get(); }
  const ChemistryContinuumOpacityTables *continuum_opacity() const noexcept { return continuum_opacity_.get(); }

private:
  Model() noexcept = default;
  ChemistryStatus initialize(const char *table_root) noexcept;
  std::unique_ptr<ChemReactionNetwork, ModelResourceDeleter> network_;
  RateTables rates_;
  std::unique_ptr<ChemistryCoolingTables, ModelResourceDeleter> molecular_cooling_;
  std::unique_ptr<ChemistryDustTables, ModelResourceDeleter> dust_opacity_;
  std::unique_ptr<ChemistryContinuumOpacityTables, ModelResourceDeleter> continuum_opacity_;
  ChemRunView view_{};
};

struct Model::Creation {
  Handle model;
  ChemistryStatus status;
};

} // namespace chemistry
