#!/bin/bash
gnuplot << 'EOF'
set terminal pngcairo size 800,600
set terminal png font 'Arial,10'
set output 'v_dist.png'
set xlabel 'v'
set ylabel 'probability'
c=0
f1(x)=a1*exp(-b1*(x-c)**2)
f2(x)=a2*exp(-b2*(x-c)**2)
f3(x)=a3*exp(-b3*(x-c)**2)
fit f1(x) 'vx_dist.dat' u 1:2 via a1,b1,c1
fit f2(x) 'vy_dist.dat' u 1:2 via a2,b2,c2
fit f3(x) 'vz_dist.dat' u 1:2 via a3,b3,c3
plot 'vx_dist.dat' u 1:2 with points linecolor rgb "blue" title 'vx', f1(x) with lines linecolor rgb "blue" title 'vx fit', 'vy_dist.dat' u 1:2 with points linecolor rgb "red" title 'vx', f2(x) with lines linecolor rgb "red" title 'vy fit', 'vz_dist.dat' u 1:2 with points linecolor rgb "green" title 'vx', f3(x) with lines linecolor rgb "green" title 'vz fit'
EOF

