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


