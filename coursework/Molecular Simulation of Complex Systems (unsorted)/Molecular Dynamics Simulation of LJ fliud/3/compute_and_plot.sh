#!/bin/bash
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

gnuplot ./plot.sh
