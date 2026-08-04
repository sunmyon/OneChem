set terminal pngcairo size 1500,1100 enhanced font "Helvetica,14"
set output "figures/rho_T_pairs_metal_J0.png"

set multiplot layout 3,2 margins 0.08,0.98,0.08,0.95 spacing 0.06,0.10
set logscale xy
set format x "10^{%L}"
set format y "10^{%L}"
set grid xtics ytics mxtics mytics
set ylabel "T [K]"

set xlabel ""
set title "Z=1e-5"
plot "figures/rho_T_comparison_Z_1em5_J0.dat" using 1:2 with lines lw 2 lc rgb "#1b9e77" title "old fixed", \
     "" using 1:3 with lines lw 2 dt 2 lc rgb "#000000" title "new"

set xlabel ""
set title "Z=1e-4"
plot "figures/rho_T_comparison_Z_1em4_J0.dat" using 1:2 with lines lw 2 lc rgb "#7570b3" title "old fixed", \
     "" using 1:3 with lines lw 2 dt 2 lc rgb "#000000" title "new"

set xlabel ""
set title "Z=1e-3"
plot "figures/rho_T_comparison_Z_1em3_J0.dat" using 1:2 with lines lw 2 lc rgb "#d95f02" title "old fixed", \
     "" using 1:3 with lines lw 2 dt 2 lc rgb "#000000" title "new"

set xlabel ""
set title "Z=1e-2"
plot "figures/rho_T_comparison_Z_1em2_J0.dat" using 1:2 with lines lw 2 lc rgb "#e7298a" title "old fixed", \
     "" using 1:3 with lines lw 2 dt 2 lc rgb "#000000" title "new"

set xlabel "n_H [cm^{-3}]"
set title "Z=1e-1"
plot "figures/rho_T_comparison_Z_1em1_J0.dat" using 1:2 with lines lw 2 lc rgb "#66a61e" title "old fixed", \
     "" using 1:3 with lines lw 2 dt 2 lc rgb "#000000" title "new"

set xlabel "n_H [cm^{-3}]"
set title "Z=1"
plot "figures/rho_T_comparison_Z_1_J0.dat" using 1:2 with lines lw 2 lc rgb "#a6761d" title "old fixed", \
     "" using 1:3 with lines lw 2 dt 2 lc rgb "#000000" title "new"

unset multiplot
