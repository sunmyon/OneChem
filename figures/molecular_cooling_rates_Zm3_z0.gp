set terminal pngcairo size 1500,900 enhanced font "Helvetica,18"
set output "figures/molecular_cooling_rates_Zm3_z0.png"
set multiplot layout 1,2 margins 0.08,0.98,0.12,0.93 spacing 0.08,0.04
set logscale xy
set xrange [1e-3:1e16]
set yrange [1e-15:1e8]
set grid xtics ytics mxtics mytics back lc rgb "#dddddd"
set key top left opaque font ",14"
set xlabel "n_{nuclei} [cm^{-3}]"
set ylabel "|du/dt| [code output units]"
set style line 1 lc rgb "#009E73" lw 2.6
set style line 2 lc rgb "#D55E00" lw 2.6
set style line 3 lc rgb "#0072B2" lw 2.6
set style line 4 lc rgb "#CC79A7" lw 2.6
set style line 5 lc rgb "#000000" lw 2.2 dt 2

set title "Molecular cooling on: Z/Z_{sun}=10^{-3}, z=0"
plot "figures/molecular_on_Zm3_z0_cool.dat" using 1:36 with lines ls 1 title "CO", \
     "figures/molecular_on_Zm3_z0_cool.dat" using 1:37 with lines ls 2 title "OH", \
     "figures/molecular_on_Zm3_z0_cool.dat" using 1:38 with lines ls 3 title "H2O", \
     "figures/molecular_on_Zm3_z0_cool.dat" using 1:28 with lines ls 4 title "dust", \
     "figures/molecular_on_Zm3_z0_cool.dat" using 1:2 with lines ls 5 title "total radiative"

set title "Molecular cooling off: same chemical network"
unset ylabel
plot "figures/molecular_off_Zm3_z0_cool.dat" using 1:36 with lines ls 1 title "CO", \
     "figures/molecular_off_Zm3_z0_cool.dat" using 1:37 with lines ls 2 title "OH", \
     "figures/molecular_off_Zm3_z0_cool.dat" using 1:38 with lines ls 3 title "H2O", \
     "figures/molecular_off_Zm3_z0_cool.dat" using 1:28 with lines ls 4 title "dust", \
     "figures/molecular_off_Zm3_z0_cool.dat" using 1:2 with lines ls 5 title "total radiative"
unset multiplot
