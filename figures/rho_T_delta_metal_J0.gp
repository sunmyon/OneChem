set terminal pngcairo size 1400,820 enhanced font "Helvetica,16"
set output "figures/rho_T_delta_metal_J0.png"

set title "1ZONE J21=0: relative gas-temperature difference after fixes"
set logscale x
set xlabel "n_H [cm^{-3}]"
set ylabel "100 x (T_new - T_old) / T_new [%]"
set grid xtics ytics mxtics
set key left bottom maxrows 2
set yrange [-0.14:0.06]
set format x "10^{%L}"

plot \
  "figures/rho_T_comparison_Z_1em5_J0.dat" using 1:($4*100.0) with lines lw 2 lc rgb "#1b9e77" title "Z=1e-5", \
  "figures/rho_T_comparison_Z_1em4_J0.dat" using 1:($4*100.0) with lines lw 2 lc rgb "#7570b3" title "Z=1e-4", \
  "figures/rho_T_comparison_Z_1em3_J0.dat" using 1:($4*100.0) with lines lw 2 lc rgb "#d95f02" title "Z=1e-3", \
  "figures/rho_T_comparison_Z_1em2_J0.dat" using 1:($4*100.0) with lines lw 2 lc rgb "#e7298a" title "Z=1e-2", \
  "figures/rho_T_comparison_Z_1em1_J0.dat" using 1:($4*100.0) with lines lw 2 lc rgb "#66a61e" title "Z=1e-1", \
  "figures/rho_T_comparison_Z_1_J0.dat" using 1:($4*100.0) with lines lw 2 lc rgb "#a6761d" title "Z=1", \
  0 with lines lw 1 lc rgb "#000000" notitle
