#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/* Lee et al. (1996) CO self-shielding and shielding by overlapping H2 lines. */
double radiation_co_line_shielding_lee1996(double column_co,
                                           double column_h2);

/* Fukushima et al. (2020), Appendix A: unshielded rate 1.03e-10 G_FUV. */
double radiation_co_photodissociation_rate(double g_fuv,
                                           double column_co,
                                           double column_h2,
                                           double dust_optical_depth);

#ifdef __cplusplus
}
#endif
