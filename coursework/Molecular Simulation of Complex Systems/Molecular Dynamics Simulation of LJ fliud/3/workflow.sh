#!/bin/bash

mkdir 2.0t 1.2t

cd ./2.0t
echo "2.0t run start"
mkdir 0.1r 0.2r 0.3r 0.4r 0.5r 0.6r 0.7r 0.8r 0.9r
for i in {1..9}; do
dir_name="0.${i}r"
cd ./${dir_name}
cp ../../compile.sh ../../md_lj.cpp ../../parameter.txt ../../prng.h ./

sed -i "3s/0.1/0.${i}/" parameter.txt
sed -i "1s/1.2/2.0/" parameter.txt

./compile.sh mdrun
./mdrun
echo "0.${i}r done!"

cd ..
done


cd ../1.2t
echo "1.2t run start"
mkdir 0.1r 0.2r 0.3r 0.4r 0.5r 0.6r 0.7r 0.8r 0.9r
for i in {1..9}; do
dir_name="0.${i}r"
cd ./${dir_name}
cp ../../compile.sh ../../md_lj.cpp ../../parameter.txt ../../prng.h ./

sed -i "3s/0.1/0.${i}/" parameter.txt
sed -i "1s/1.2/1.2/" parameter.txt
sed -i "4s/0.001/0.0001/" parameter.txt

./compile.sh mdrun
./mdrun
echo "0.${i}r done!"

cd ..
done




cd ./1.2t
touch all_p.txt

for i in {1..9}; do
dir_name="0.${i}r"

awk -v i="$i" '{print "0."i, $0}' ./${dir_name}/pressure.dat >> all_p.txt


done

cd ../2.0t
touch all_p.txt

for i in {1..9}; do
dir_name="0.${i}r"

awk -v i="$i" '{print "0."i, $0}' ./${dir_name}/pressure.dat >> all_p.txt


done

cd ..




cd ./1.2t

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
}' all_p.txt > comp_p.txt

sort -n -k1,1 -o comp_p.txt comp_p.txt

cd ../2.0t
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
}' all_p.txt > comp_p.txt

sort -n -k1,1 -o comp_p.txt comp_p.txt

cd ..



gnuplot << 'EOF'

set terminal pngcairo enhanced
set output 'p_T.png'

set title "MD simulation: p V.S. T"
set xlabel "T"
set ylabel "pressure"
set grid

set yrange [*:*]

plot "./1.2t/comp_p.txt" using 1:2:3 with yerrorlines title "T=1.2" lt 1 lc rgb "blue" pt 7 lw 2, \
"./2.0t/comp_p.txt" using 1:2:3 with yerrorlines title "T=2.0" lt 1 lc rgb "red" pt 7 lw 2



EOF
