set terminal pngcairo size 1400,980 enhanced font "Helvetica,15"
set output "figures/rho_T_comparison_metal_J0.png"

set multiplot layout 2,1 margins 0.10,0.97,0.10,0.94 spacing 0.03,0.02

set title "1ZONE J21=0: metallicity sweep rho-T comparison (old fixed, adapter dust tau fixed)"
set logscale xy
set format x ""
set ylabel "Temperature [K]"
set grid xtics ytics mxtics mytics
set key left top maxrows 6
plot \
  "figures/rho_T_comparison_Z_1em5_J0.dat" using 1:2 with lines lw 2 lc rgb "#1b9e77" title "old fixed Z=1e-5", \
  "" using 1:3 with lines lw 2 dt 2 lc rgb "#1b9e77" title "new Z=1e-5", \
  "figures/rho_T_comparison_Z_1em4_J0.dat" using 1:2 with lines lw 2 lc rgb "#7570b3" title "old fixed Z=1e-4", \
  "" using 1:3 with lines lw 2 dt 2 lc rgb "#7570b3" title "new Z=1e-4", \
  "figures/rho_T_comparison_Z_1em3_J0.dat" using 1:2 with lines lw 2 lc rgb "#d95f02" title "old fixed Z=1e-3", \
  "" using 1:3 with lines lw 2 dt 2 lc rgb "#d95f02" title "new Z=1e-3", \
  "figures/rho_T_comparison_Z_1em2_J0.dat" using 1:2 with lines lw 2 lc rgb "#e7298a" title "old fixed Z=1e-2", \
  "" using 1:3 with lines lw 2 dt 2 lc rgb "#e7298a" title "new Z=1e-2", \
  "figures/rho_T_comparison_Z_1em1_J0.dat" using 1:2 with lines lw 2 lc rgb "#66a61e" title "old fixed Z=1e-1", \
  "" using 1:3 with lines lw 2 dt 2 lc rgb "#66a61e" title "new Z=1e-1", \
  "figures/rho_T_comparison_Z_1_J0.dat" using 1:2 with lines lw 2 lc rgb "#a6761d" title "old fixed Z=1", \
  "" using 1:3 with lines lw 2 dt 2 lc rgb "#a6761d" title "new Z=1"

unset title
unset logscale y
set logscale x
set format x "10^{%L}"
set xlabel "n_H [cm^{-3}]"
set ylabel "{/Symbol D}T/T_{new}"
set grid xtics ytics mxtics
set key left top maxrows 3
plot \
  "figures/rho_T_comparison_Z_1em5_J0.dat" using 1:4 with lines lw 2 lc rgb "#1b9e77" title "1e-5", \
  "figures/rho_T_comparison_Z_1em4_J0.dat" using 1:4 with lines lw 2 lc rgb "#7570b3" title "1e-4", \
  "figures/rho_T_comparison_Z_1em3_J0.dat" using 1:4 with lines lw 2 lc rgb "#d95f02" title "1e-3", \
  "figures/rho_T_comparison_Z_1em2_J0.dat" using 1:4 with lines lw 2 lc rgb "#e7298a" title "1e-2", \
  "figures/rho_T_comparison_Z_1em1_J0.dat" using 1:4 with lines lw 2 lc rgb "#66a61e" title "1e-1", \
  "figures/rho_T_comparison_Z_1_J0.dat" using 1:4 with lines lw 2 lc rgb "#a6761d" title "1", \
  0 with lines lw 1 lc rgb "#000000" notitle

unset multiplot
