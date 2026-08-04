set terminal pngcairo size 1200,900 enhanced font "Helvetica,16"
set output "figures/rho_T_comparison_nolw_0metal.png"

set multiplot layout 2,1 margins 0.10,0.97,0.10,0.94 spacing 0.03,0.02

set title "1ZONE no-LW, 0-metal: rho-T comparison"
set logscale xy
set format x ""
set ylabel "Temperature [K]"
set grid xtics ytics mxtics mytics
set key left top
plot "figures/rho_T_comparison_nolw_0metal.dat" using 1:2 with lines lw 3 lc rgb "#2c5aa0" title "chemistry_old", \
     "" using 1:3 with lines lw 3 dt 2 lc rgb "#d95f02" title "chemistry"

unset title
unset logscale y
set logscale x
set format x "10^{%L}"
set xlabel "n_H [cm^{-3}]"
set ylabel "{/Symbol D}T/T_{old}"
set grid xtics ytics mxtics
unset key
plot "figures/rho_T_comparison_nolw_0metal.dat" using 1:4 with lines lw 2 lc rgb "#444444", \
     0 with lines lw 1 lc rgb "#000000"

unset multiplot
