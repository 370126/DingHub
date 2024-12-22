#!/bin/bash
gnuplot << 'EOF'
fileinp='vx_dist.dat'
set terminal png font 'Arial,10'
set output 'vx_dist.png'
set xlabel 'vx'
set ylabel 'probability'
c=0
f(x)=a*exp(-b*(x-c)**2)
fit f(x) fileinp u 1:2 via a,b,c
plot fileinp u 1:2 w p pt 6 ps 2 lc rgb "black" not, f(x) not lt 7 lc rgb "blue"
EOF

gnuplot << 'EOF'
fileinp='vy_dist.dat'
set terminal png font 'Arial,10'
set output 'vy_dist.png'
set xlabel 'vy'
set ylabel 'probability'
c=0
f(x)=a*exp(-b*(x-c)**2)
fit f(x) fileinp u 1:2 via a,b,c
plot fileinp u 1:2 w p pt 6 ps 2 lc rgb "black" not, f(x) not lt 7 lc rgb "blue"
EOF

gnuplot << 'EOF'
fileinp='vz_dist.dat'
set terminal png font 'Arial,10'
set output 'vz_dist.png'
set xlabel 'vz'
set ylabel 'probability'
c=1
f(x)=a*exp(-b*(x-c)**2)
fit f(x) fileinp u 1:2 via a,b,c
plot fileinp u 1:2 w p pt 6 ps 2 lc rgb "black" not, f(x) not lt 7 lc rgb "blue"
EOF

gnuplot << 'EOF'
fileinp='v_dist.dat'
set terminal png font 'Arial,10'
set output 'v_dist.png'
set xlabel 'v'
set ylabel 'probability'
c=2
f(x)=a*exp(-b*(x-c)**2)
fit f(x) fileinp u 1:2 via a,b,c
plot fileinp u 1:2 w p pt 6 ps 2 lc rgb "black" not, f(x) not lt 7 lc rgb "blue"
EOF
