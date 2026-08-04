set terminal pngcairo size 1200,520 enhanced font "Helvetica,14"
set output "figures/metal_lw_heatmap_g0unit.png"
set multiplot layout 1,2 margins 0.08,0.96,0.14,0.88 spacing 0.10,0.02
set pm3d map
set view map
set datafile separator whitespace
set xlabel "log10 J21"
set ylabel "log10 Z/Zsun"
set xrange [-1.2:3.2]
set yrange [-5.3:0.3]
set xtics ("0.1" -1, "1" 0, "10" 1, "100" 2, "1000" 3)
set ytics ("1e-5" -5, "1e-4" -4, "1e-3" -3, "1e-2" -2, "1e-1" -1, "1" 0)
set cblabel "relative difference"
set logscale cb
set cbrange [1e-8:1e-3]
set palette defined (0 "#17324d", 0.35 "#2b8cbe", 0.65 "#fdae61", 1 "#9e0142")
set title "max |dT/T| over trajectory"
splot "figures/metal_lw_heatmap_g0unit.dat" using 1:2:3 with pm3d notitle
set title "final |dT/T|"
splot "figures/metal_lw_heatmap_g0unit.dat" using 1:2:(abs($4)) with pm3d notitle
unset multiplot
