#!/bin/bash

rm all_p1_ave_rho.txt
touch all_p1_ave_rho.txt

#extract
for f in {1..7}; do

dir=$f

cd $dir/ 

#./compile.sh
#./dpd

awk -v i="$f" '{print (($1+$5+$9)/3-i)/25}' pressure_tensor.out > p1_ave.txt
awk -v i="$f" '{print i, $0}' p1_ave.txt > p1_ave_rho.txt 

echo "done in $f"

cat p1_ave_rho.txt >> ./../all_p1_ave_rho.txt

cd ..

done

# compute
awk '{
x[$1] += $2
x_sq[$1] += $2^2
count[$1]++
}
END {
for (i in x) {
	mean = x[i] / count[i]
	va = (x_sq[i] / count[i]) - (mean^2)
	sd = sqrt(va)
	printf "%s %f %f\n", i, mean, sd
	}
}' all_p1_ave_rho.txt > comp_p1.txt

sort -n -k1,1 -o comp_p1.txt comp_p1.txt


# plot
gnuplot << 'EOF'

set terminal pngcairo enhanced size 1000,600
set output 'p1_rho.png'

set title "DPD simulation: (p-rho*kT)/a V.S. rho"
set xlabel "rho"
set ylabel "(p-rho*kT)/a"
set grid

set yrange [*:*]
set xrange [*:*]

plot "comp_p1.txt" using 1:2:3 with yerrorlines title "p average | std" lt 1 pt 5 lw 2 lc rgb "black"



EOF
