awk '{print $2,$3,$4}' traj.xyz > xyz.dat
awk '{print $1+$2}' ener_momentum.dat > tot_ener.dat
awk '{print NR, $1}' tot_ener.dat > plot_ener.txt
awk '{if((NR-1)%502>1){print $5,$6,$7,sqrt($5**2+$6**2+$7**2)}}' traj.xyz > vel.dat


./plot_energy



