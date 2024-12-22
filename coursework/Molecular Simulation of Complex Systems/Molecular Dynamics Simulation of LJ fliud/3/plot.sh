

set terminal pngcairo enhanced
set output 'p_T.png'

set title "MD simulation: p V.S. T"
set xlabel "T"
set ylabel "pressure"
set grid

set yrange [*:*]

plot "./1.2t/comp_p.txt" using 1:2:3 with yerrorlines title "T=1.2" lt 1 lc rgb "blue" pt 7 lw 2, \
"./2.0t/comp_p.txt" using 1:2:3 with yerrorlines title "T=2.0" lt 1 lc rgb "red" pt 7 lw 2

