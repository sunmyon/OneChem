set terminal pngcairo size 1300,950 enhanced font "Helvetica,15"
set output "figures/rho_T_comparison_lw_0metal.png"

set multiplot layout 2,1 margins 0.10,0.97,0.10,0.94 spacing 0.03,0.02

set title "1ZONE 0-metal with background LW: rho-T comparison"
set logscale xy
set format x ""
set ylabel "Temperature [K]"
set grid xtics ytics mxtics mytics
set key left top maxrows 5
plot \
  "figures/rho_T_comparison_J21_0p1_0metal.dat" using 1:2 with lines lw 2 lc rgb "#1b9e77" title "old J21=0.1", \
  "" using 1:3 with lines lw 2 dt 2 lc rgb "#1b9e77" title "new J21=0.1", \
  "figures/rho_T_comparison_J21_1_0metal.dat" using 1:2 with lines lw 2 lc rgb "#7570b3" title "old J21=1", \
  "" using 1:3 with lines lw 2 dt 2 lc rgb "#7570b3" title "new J21=1", \
  "figures/rho_T_comparison_J21_10_0metal.dat" using 1:2 with lines lw 2 lc rgb "#d95f02" title "old J21=10", \
  "" using 1:3 with lines lw 2 dt 2 lc rgb "#d95f02" title "new J21=10", \
  "figures/rho_T_comparison_J21_100_0metal.dat" using 1:2 with lines lw 2 lc rgb "#e7298a" title "old J21=100", \
  "" using 1:3 with lines lw 2 dt 2 lc rgb "#e7298a" title "new J21=100", \
  "figures/rho_T_comparison_J21_1000_0metal.dat" using 1:2 with lines lw 2 lc rgb "#66a61e" title "old J21=1000", \
  "" using 1:3 with lines lw 2 dt 2 lc rgb "#66a61e" title "new J21=1000"

unset title
unset logscale y
set logscale x
set format x "10^{%L}"
set xlabel "n_H [cm^{-3}]"
set ylabel "{/Symbol D}T/T_{old}"
set grid xtics ytics mxtics
set key left top maxrows 3
plot \
  "figures/rho_T_comparison_J21_0p1_0metal.dat" using 1:4 with lines lw 2 lc rgb "#1b9e77" title "0.1", \
  "figures/rho_T_comparison_J21_1_0metal.dat" using 1:4 with lines lw 2 lc rgb "#7570b3" title "1", \
  "figures/rho_T_comparison_J21_10_0metal.dat" using 1:4 with lines lw 2 lc rgb "#d95f02" title "10", \
  "figures/rho_T_comparison_J21_100_0metal.dat" using 1:4 with lines lw 2 lc rgb "#e7298a" title "100", \
  "figures/rho_T_comparison_J21_1000_0metal.dat" using 1:4 with lines lw 2 lc rgb "#66a61e" title "1000", \
  0 with lines lw 1 lc rgb "#000000" notitle

unset multiplot
