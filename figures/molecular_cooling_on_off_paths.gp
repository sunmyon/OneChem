set terminal pngcairo size 1500,1050 enhanced font "Helvetica,18"
set output "figures/molecular_cooling_on_off_paths.png"
set multiplot layout 2,2 rowsfirst margins 0.09,0.98,0.09,0.95 spacing 0.06,0.12
set logscale xy
set xrange [1e-3:1e16]
set yrange [3:1e4]
set grid xtics ytics mxtics mytics back lc rgb "#dddddd"
set key top right opaque
set ylabel "T_{gas} [K]"
set style line 1 lc rgb "#0072B2" lw 3
set style line 2 lc rgb "#D55E00" lw 3 dt 2
set style line 3 lc rgb "#009E73" lw 3 dt 3

set title "Z/Z_{sun} = 10^{-3}, z = 0"
unset xlabel
plot "figures/molecular_on_Zm3_z0.dat" using 1:2 with lines ls 1 title "molecular cooling on", \
     "figures/molecular_off_Zm3_z0.dat" using 1:2 with lines ls 2 title "molecular cooling off", \
     "figures/no_co_chem_CII_OI_Zm3_z0.dat" using 1:2 with lines ls 3 title "C/O chemistry off; CII+OI cooling"

set title "Z/Z_{sun} = 10^{-2}, z = 0"
unset ylabel
plot "figures/molecular_on_Zm2_z0.dat" using 1:2 with lines ls 1 title "on", \
     "figures/molecular_off_Zm2_z0.dat" using 1:2 with lines ls 2 title "mol. cooling off", \
     "figures/no_co_chem_CII_OI_Zm2_z0.dat" using 1:2 with lines ls 3 title "C/O chemistry off"

set title "Z/Z_{sun} = 10^{-1}, z = 0"
set xlabel "n_{nuclei} [cm^{-3}]"
set ylabel "T_{gas} [K]"
plot "figures/molecular_on_Zm1_z0.dat" using 1:2 with lines ls 1 title "on", \
     "figures/molecular_off_Zm1_z0.dat" using 1:2 with lines ls 2 title "mol. cooling off", \
     "figures/no_co_chem_CII_OI_Zm1_z0.dat" using 1:2 with lines ls 3 title "C/O chemistry off"

set title "Z/Z_{sun} = 1, z = 0"
unset ylabel
plot "figures/molecular_on_Z0_z0.dat" using 1:2 with lines ls 1 title "on", \
     "figures/molecular_off_Z0_z0.dat" using 1:2 with lines ls 2 title "mol. cooling off", \
     "figures/no_co_chem_CII_OI_Z0_z0.dat" using 1:2 with lines ls 3 title "C/O chemistry off"
unset multiplot
