#!/bin/bash
awk '{print $1,$2/2/($1+0.0000000000001)}' diff.dat > msd_x.dat
awk '{print $1,$3/2/($1+0.0000000000001)}' diff.dat > msd_y.dat
awk '{print $1,$4/2/($1+0.0000000000001)}' diff.dat > msd_z.dat

gnuplot << 'EOF'

set terminal pngcairo size 1000,600 enhanced font 'Arial,10'
set output 'DEL_MSD_vs_t.png'
set xlabel 't'
set ylabel '⟨∆rx^2(t)⟩/2dt'
set yrange [*:*]
set grid

plot 'msd_x.dat' using 1:2 with lines title 'Dx' lc rgb 'blue', \
'msd_y.dat' using 1:2 with lines title 'Dy' lc rgb 'red', \
'msd_z.dat' using 1:2 with lines title 'Dz' lc rgb 'green'
EOF



# GOOD
awk '$1>45 && $1<55' msd_x.dat > good_msd_x.dat
awk '$1>45 && $1<55' msd_y.dat > good_msd_y.dat
awk '$1>45 && $1<55' msd_z.dat > good_msd_z.dat

awk '{sum+=$2} END {print sum/NR}' good_msd_x.dat > D_x.dat
awk '{sum+=$2} END {print sum/NR}' good_msd_y.dat > D_y.dat
awk '{sum+=$2} END {print sum/NR}' good_msd_z.dat > D_z.dat


cat good_msd_*.dat | awk '{sum+=$2} END {print sum/NR}' > D_ave_MSD.dat


gnuplot << 'EOF'

set terminal pngcairo size 1000,600 enhanced font 'Arial,10'
set output 'GOOD_DEL_MSD_vs_t.png'
set xlabel 't'
set ylabel '⟨∆rx^2(t)⟩/2d'
set yrange [*:*]
set xrange [*:*]
set grid

plot 'good_msd_x.dat' using 1:2 with lines title 'Dx' lc rgb 'blue', \
'good_msd_y.dat' using 1:2 with lines title 'Dy' lc rgb 'red', \
'good_msd_z.dat' using 1:2 with lines title 'Dz' lc rgb 'green'
EOF


rm good_msd_*.dat
rm msd_*.dat
rm D_x.dat
rm D_y.dat
rm D_z.dat
