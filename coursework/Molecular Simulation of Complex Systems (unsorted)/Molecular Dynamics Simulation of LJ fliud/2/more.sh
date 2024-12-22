awk '{print $2,$3,$4}' traj.xyz > xyz.dat
awk '{print $1+$2}' ener_momentum.dat > tot_ener.dat
awk '{print NR, $1}' tot_ener.dat > plot_ener_eular.txt
awk '{print $3}' ener_momentum.dat > vx.dat
awk '{print $4}' ener_momentum.dat > vy.dat
awk '{print $5}' ener_momentum.dat > vz.dat
awk '{print sqrt($3*$3+$4*$4+$5*$5)}' ener_momentum.dat > v.dat



gnuplot plot



