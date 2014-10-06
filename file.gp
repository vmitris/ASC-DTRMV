set term png
set output "opteron.png"
set title "Opteron"
set xlabel "Dimensiune"
set ylabel "Timp(s)"
set grid
set style line 22 linetype 3 linewidth 4
plot "opteron.out" using 1:2 title "Rulare de mana" with linespoints, "opteron.out" using 1:3 title "Rulare Blas" with linespoints, "opteron.out" using 1:4 title "Optimizat" with linespoints
set output "nehalem.png"
set title "Nehalem"
plot "nehalem.out" using 1:2 title "Rulare de mana" with linespoints, "nehalem.out" using 1:3 title "Rulare Blas" with linespoints, "nehalem.out" using 1:4 title "Optimizat" with linespoints
set output "quad.png"
set title "Quad"
plot "quad.out" using 1:2 title "Rulare de mana" with linespoints, "quad.out" using 1:3 title "Rulare Blas" with linespoints, "quad.out" using 1:4 title "Optimizat" with linespoints



