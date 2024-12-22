#!/bin/bash

rm -f D_acf_xyz.dat
touch D_acf_xyz.dat
awk '$1<50 {sum+=$5}END{print sum*0.01}' diff.dat >> D_acf_xyz.dat
awk '$1<50 {sum+=$6}END{print sum*0.01}' diff.dat >> D_acf_xyz.dat
awk '$1<50 {sum+=$7}END{print sum*0.01}' diff.dat >> D_acf_xyz.dat

awk '{sum+=$1}END{print sum/NR}' D_acf_xyz.dat > D_ave_ACF.dat



gnuplot << 'EOF'

set terminal pngcairo size 1000,600 enhanced 
set output 'ACF vs t.png'
set xlabel 't'
set ylabel 'Integrated ACF'
set yrange [*:*]
set grid

plot 'diff.dat' using 1:5 with lines title 'ACF_x' lc rgb 'blue', \
'diff.dat' using 1:6 with lines title 'ACF_y' lc rgb 'red', \
'diff.dat' using 1:7 with lines title 'ACF_z' lc rgb 'green'
EOF
