#!/bin/bash
#g++ -o2 DPD.cpp -o $1 -DTENSION
g++ -o2 -fopenmp DPD.cpp -pthread -o $1 -DTENSION
