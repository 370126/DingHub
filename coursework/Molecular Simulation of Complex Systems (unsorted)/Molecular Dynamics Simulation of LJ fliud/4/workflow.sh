#!/bin/bash

mkdir 0.1r 0.2r 0.3r 0.4r 0.5r 0.6r 0.7r 0.8r 0.9r
for i in {1..9}; do
dir_name="0.${i}r"
cd ./${dir_name}
cp ../*.sh ../*.cpp  ../*.txt ../*.h ./

sed -i "3s/0.1/0.${i}/" parameter.txt
echo "0.${i}r run start"
./compile.sh mdrun
./mdrun
echo "0.${i}r MD run done!"

./compile_diff.sh diff_run
./diff_run
echo "0.${i}r Diff run done!"

./MSD.sh
echo "MSD calculation done"
./ACF.sh
echo "ACF calculation done"

cd ..
echo "0.${i}r run done!"
done


