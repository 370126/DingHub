#!/bin/bash

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
