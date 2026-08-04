set terminal pngcairo size 1600,760 enhanced font "Helvetica,17"
set output "figures/three_model_performance_comparison.png"
set multiplot layout 1,2 margins 0.07,0.98,0.13,0.91 spacing 0.09,0.04
set style data histograms
set style histogram rowstacked
set style fill solid 0.85 border -1
set boxwidth 0.75
set grid ytics back lc rgb "#dddddd"
set yrange [0:114]
set ylabel "fraction of chemistry wall time [%]"
set key inside top center horizontal font ",13"
set title "Exclusive top-level time breakdown"
plot "figures/three_model_performance_comparison.dat" using 2:xtic(1) title "implicit chemistry" lc rgb "#0072B2", \
     "" using 3 title "radiative cooling" lc rgb "#D55E00", \
     "" using 4 title "residual overhead" lc rgb "#009E73", \
     "" using 5 title "solver + driver" lc rgb "#CC79A7"

set style histogram clustered gap 1
set yrange [0:6.8]
set ylabel "mean evaluations / parent call"
set key inside top center horizontal font ",13"
set title "Iterations needed for convergence"
plot "figures/three_model_performance_comparison.dat" using 6:xtic(1) title "energy trials / collapse step" lc rgb "#56B4E9", \
     "" using 7 title "Newton iterations / implicit trial (all)" lc rgb "#E69F00"
unset multiplot
