#!/bin/bash

rm -f data_all.txt
touch data_all.txt
for i in {1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.7}; do

awk -v i="$i" '{print i, $0}' ./${i}/Sigma_result.txt >> data_all.txt

done



# plot
gnuplot << 'EOF'

set terminal pngcairo enhanced size 1000,650
set output 'Sigma_A.png'

set title "DPD simulation: "
set xlabel "A"
set ylabel "Sigma"
set grid

set yrange [*:*]
set xrange [1.1:1.8]

plot "data_all.txt" using 1:2 with linespoints title "Sigma\_ave" ps 2 pt 1 lw 1.25 lc rgb "blue"

EOF
