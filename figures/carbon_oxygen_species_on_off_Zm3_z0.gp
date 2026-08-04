set terminal pngcairo size 1550,1100 enhanced font "Helvetica,17"
set output "figures/carbon_oxygen_species_on_off_Zm3_z0.png"
set multiplot layout 2,2 rowsfirst margins 0.09,0.98,0.09,0.95 spacing 0.06,0.12
set logscale xy
set xrange [1e-3:1e16]
set yrange [1e-20:1e-6]
set grid xtics ytics mxtics mytics back lc rgb "#dddddd"
set key bottom left opaque font ",14"
set ylabel "abundance / X_H"
set style line 1 lc rgb "#000000" lw 2.6
set style line 2 lc rgb "#E69F00" lw 2.6
set style line 3 lc rgb "#009E73" lw 2.6
set style line 4 lc rgb "#CC79A7" lw 2.6

set title "Carbon species: molecular cooling on"
unset xlabel
plot "figures/molecular_on_Zm3_z0.dat" using 1:30 with lines ls 1 title "C", \
     "figures/molecular_on_Zm3_z0.dat" using 1:31 with lines ls 2 title "C+", \
     "figures/molecular_on_Zm3_z0.dat" using 1:32 with lines ls 3 title "CO", \
     "figures/molecular_on_Zm3_z0.dat" using 1:33 with lines ls 4 title "CO2"

set title "Carbon species: molecular cooling off"
unset ylabel
plot "figures/molecular_off_Zm3_z0.dat" using 1:30 with lines ls 1 title "C", \
     "figures/molecular_off_Zm3_z0.dat" using 1:31 with lines ls 2 title "C+", \
     "figures/molecular_off_Zm3_z0.dat" using 1:32 with lines ls 3 title "CO", \
     "figures/molecular_off_Zm3_z0.dat" using 1:33 with lines ls 4 title "CO2"

set style line 5 lc rgb "#0072B2" lw 2.6
set style line 6 lc rgb "#D55E00" lw 2.6
set style line 7 lc rgb "#56B4E9" lw 2.6
set style line 8 lc rgb "#F0E442" lw 2.6
set title "Oxygen species: molecular cooling on"
set xlabel "n_{nuclei} [cm^{-3}]"
set ylabel "abundance / X_H"
plot "figures/molecular_on_Zm3_z0.dat" using 1:34 with lines ls 5 title "O", \
     "figures/molecular_on_Zm3_z0.dat" using 1:35 with lines ls 6 title "OH", \
     "figures/molecular_on_Zm3_z0.dat" using 1:36 with lines ls 7 title "H2O", \
     "figures/molecular_on_Zm3_z0.dat" using 1:37 with lines ls 8 title "O2"

set title "Oxygen species: molecular cooling off"
unset ylabel
plot "figures/molecular_off_Zm3_z0.dat" using 1:34 with lines ls 5 title "O", \
     "figures/molecular_off_Zm3_z0.dat" using 1:35 with lines ls 6 title "OH", \
     "figures/molecular_off_Zm3_z0.dat" using 1:36 with lines ls 7 title "H2O", \
     "figures/molecular_off_Zm3_z0.dat" using 1:37 with lines ls 8 title "O2"
unset multiplot
