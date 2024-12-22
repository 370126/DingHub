#!/bin/bash

rm -f D_MSD_all.dat
rm -f D_ACF_all.dat
touch D_MSD_all.dat
touch D_ACF_all.dat

for i in {1..9}; do
dir="0.${i}r"
rho="0.${i}"
D_msd=$(cat $dir/D_ave_MSD.dat)
D_acf=$(cat $dir/D_ave_ACF.dat)
echo "$rho	$D_msd" >> D_MSD_all.dat
echo "$rho	$D_acf" >> D_ACF_all.dat
done

gnuplot << 'EOF'

set terminal pngcairo size 1000,600 enhanced 
set output 'D_all.png'

set title "MD simulation: D v.s. rho"
set xlabel "rho"
set ylabel "D"
set grid

set yrange [*:*]

plot "./D_MSD_all.dat" using 1:2 with linespoints title "from Einstein relation" lt 1 lc rgb "blue" pt 8 lw 2, \
"./D_ACF_all.dat" using 1:2 with linespoints title "Green-Kubo relation" lt 1 lc rgb "red" pt 4 lw 2

EOF
