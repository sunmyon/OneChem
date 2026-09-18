#include "chemistry_models.h"

#if defined(H2_COOLING_G15) + defined(H2_COOLING_GP12) + defined(H2_COOLING_Y06) + defined(H2_COOLING_O00) > 1
#error "Select only one H2 cooling model"
#endif

#if defined(PARTITION_H2_EBM1989) + defined(PARTITION_H2_I87) + defined(PARTITION_H2_R85) > 1
#error "Select only one H2 partition-function model"
#endif

#if defined(H2_COOLING_Y06) || defined(H2_COOLING_O00)
#error "The selected H2 cooling model has no implementation in the portable chemistry core"
#elif defined(H2_COOLING_GP12)
#define SELECTED_H2_COOLING H2_cooling_GP12
#else
#define SELECTED_H2_COOLING H2_cooling_G15
#endif

#ifdef PARTITION_H2_EBM1989
#define SELECTED_H2_PARTITION partition_H2_EBM1989
#elif defined(PARTITION_H2_I87)
#define SELECTED_H2_PARTITION partition_H2_Irwin1987
#elif defined(PARTITION_H2_R85)
#define SELECTED_H2_PARTITION partition_H2_Rossi1985
#else
#define SELECTED_H2_PARTITION partition_H2_EBM1989_fitted_by_Matsukoba
#endif

#ifdef H2_OPACITY
#ifdef H2_COOLING_O00
#define SELECTED_H2_LINE_ESCAPE fesc_unity
#elif defined(H2_OPACITY_FIT)
#define SELECTED_H2_LINE_ESCAPE fesc_fit
#else
#error "H2_OPACITY requires an implemented H2 line-escape model"
#endif
#ifdef HD
#ifdef H2_COOLING_O00
#define SELECTED_HD_LINE_ESCAPE fesc_unity
#else
#define SELECTED_HD_LINE_ESCAPE fesc_HD_fit
#endif
#endif
#endif

const ChemistryModels chemistry_models = {
  .h2_cooling = SELECTED_H2_COOLING,
  .h2_partition = SELECTED_H2_PARTITION,
#ifdef H2_OPACITY
  .h2_line_escape = SELECTED_H2_LINE_ESCAPE,
#ifdef HD
  .hd_line_escape = SELECTED_HD_LINE_ESCAPE,
#endif
#endif
};
