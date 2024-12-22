#!/bin/bash
awk 'BEGIN{NBINS=20;ndata=0;max=-10000;min=10000;}\
     {ndata++;ener[ndata]=$1;if(max<$1){max=$1};if(min>$1){min=$1}}\
     END{bin=(max-min)/NBINS;for(i=1;i<=NR;i++){n_hist=int((ener[i]-min)/bin)+1;if(n_hist>NBINS){n_hist=NBINS};count[n_hist]++};\
     for(n=1;n<=NBINS;n++){print min+bin*(n-0.5),count[n]/NR}}' vel.dat > vx_dist.dat

awk 'BEGIN{NBINS=20;ndata=0;max=-10000;min=10000;}\
     {ndata++;ener[ndata]=$2;if(max<$2){max=$2};if(min>$2){min=$2}}\
     END{bin=(max-min)/NBINS;for(i=1;i<=NR;i++){n_hist=int((ener[i]-min)/bin)+1;if(n_hist>NBINS){n_hist=NBINS};count[n_hist]++};\
     for(n=1;n<=NBINS;n++){print min+bin*(n-0.5),count[n]/NR}}' vel.dat > vy_dist.dat

awk 'BEGIN{NBINS=20;ndata=0;max=-10000;min=10000;}\
     {ndata++;ener[ndata]=$3;if(max<$3){max=$3};if(min>$3){min=$3}}\
     END{bin=(max-min)/NBINS;for(i=1;i<=NR;i++){n_hist=int((ener[i]-min)/bin)+1;if(n_hist>NBINS){n_hist=NBINS};count[n_hist]++};\
     for(n=1;n<=NBINS;n++){print min+bin*(n-0.5),count[n]/NR}}' vel.dat > vz_dist.dat

awk 'BEGIN{NBINS=20;ndata=0;max=-10000;min=10000;}\
     {ndata++;ener[ndata]=$4;if(max<$4){max=$4};if(min>$4){min=$4}}\
     END{bin=(max-min)/NBINS;for(i=1;i<=NR;i++){n_hist=int((ener[i]-min)/bin)+1;if(n_hist>NBINS){n_hist=NBINS};count[n_hist]++};\
     for(n=1;n<=NBINS;n++){print min+bin*(n-0.5),count[n]/NR}}' vel.dat > v_dist.dat
