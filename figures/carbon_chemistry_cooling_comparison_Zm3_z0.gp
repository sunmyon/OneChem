set terminal pngcairo size 1800,650 enhanced font "Helvetica,17"
set output "figures/carbon_chemistry_cooling_comparison_Zm3_z0.png"
set multiplot layout 1,3 margins 0.07,0.99,0.14,0.91 spacing 0.05,0.03
set logscale xy
set xrange [1e-3:1e16]
set yrange [1e-15:1e8]
set grid xtics ytics mxtics mytics back lc rgb "#dddddd"
set xlabel "n_{nuclei} [cm^{-3}]"
set ylabel "|du/dt| [code output units]"
set key top left opaque font ",13"
set style line 1 lc rgb "#E69F00" lw 2.5
set style line 2 lc rgb "#0072B2" lw 2.5
set style line 3 lc rgb "#009E73" lw 2.5
set style line 4 lc rgb "#D55E00" lw 2.5
set style line 5 lc rgb "#CC79A7" lw 2.5
set style line 6 lc rgb "#000000" lw 2.2 dt 2

set title "C/O chemistry + molecular cooling"
plot "figures/molecular_on_Zm3_z0_cool.dat" using 1:26 with lines ls 1 title "CII", \
     "figures/molecular_on_Zm3_z0_cool.dat" using 1:27 with lines ls 2 title "OI", \
     "figures/molecular_on_Zm3_z0_cool.dat" using 1:36 with lines ls 3 title "CO", \
     "figures/molecular_on_Zm3_z0_cool.dat" using 1:37 with lines ls 4 title "OH", \
     "figures/molecular_on_Zm3_z0_cool.dat" using 1:38 with lines ls 5 title "H2O", \
     "figures/molecular_on_Zm3_z0_cool.dat" using 1:2 with lines ls 6 title "total"

unset ylabel
set title "C/O chemistry; molecular cooling off"
plot "figures/molecular_off_Zm3_z0_cool.dat" using 1:26 with lines ls 1 title "CII", \
     "figures/molecular_off_Zm3_z0_cool.dat" using 1:27 with lines ls 2 title "OI", \
     "figures/molecular_off_Zm3_z0_cool.dat" using 1:2 with lines ls 6 title "total"

set title "C/O chemistry off; fixed CII + OI"
plot "figures/no_co_chem_CII_OI_Zm3_z0_cool.dat" using 1:26 with lines ls 1 title "CII", \
     "figures/no_co_chem_CII_OI_Zm3_z0_cool.dat" using 1:27 with lines ls 2 title "OI", \
     "figures/no_co_chem_CII_OI_Zm3_z0_cool.dat" using 1:2 with lines ls 6 title "total"
unset multiplot
