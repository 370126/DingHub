#!/bin/bash

#mkdir 1.3 1.35 1.4 1.45 1.5 1.55 1.6
mkdir 1.7
for i in {1.25,1.7}; do

cd ${i}/
cp ../*.h ../compile.sh ../*.cpp ./
sed -i "67s/1.*/${i}/" DPD.cpp
echo "$i run start"
./compile.sh dpdrun
./dpdrun
awk '{av+=$10}END{print av/NR}' stress_tensor.out > Sigma_result.txt
cd ..
done
