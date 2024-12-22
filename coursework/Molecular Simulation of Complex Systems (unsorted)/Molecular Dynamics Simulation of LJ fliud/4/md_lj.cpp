/**********************************************************
 NVE Molecular Dynamics simulations of Lennard-Jones fluid
 written by HU, Yi & HU, Jinglei 2016.03.23
**********************************************************/
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <cmath>
#include "prng.h"

#define PI 3.141592653589793238464

/* global variables in the simulations */
cPRNG<KISS, UNIFORM> prng; // random number generator
int    n_particles; // total number of LJ particles
double temp, beta, mass = 1.0; // beta = 1.0 / kT
double rho, vol, box, half_box, lbox[3], half_lbox[3]; // density, volume and simulation box size
double lj_rcut, lj_rcut_sq, lj_sigma, lj_epsilon, lj_ucut, lj_utail, lj_ptail; // LJ parameters
double **pos, **vel, **fff; // 2D array for particle positions (xi,yi,zi), (vxi, vyi, vzi), (fxi, fyi, fzi)
int    **pbc;
int    relax_steps, prod_steps, freq_ener = 100, freq_sample, freq_scale; // relaxation and equilibrium steps
double md_delt, md_deltsq, md_delt_half;      // delta t of integrate

/* Verlet lists */
int    max_num_neigh, *num_neigh, **neigh_list;
int    nl_update_count = 0;    // maximum verlet lists length, verlet list update count
double nl_rcut, nl_rcut_sq, nl_rskin, **dis;         // verlet list sphere radii, and the shell thickness

/* functions */
inline void apply_min_img(double* r) { // minimum image convention
  for (int dim = 0; dim < 3; ++dim) {
    if (r[dim] > half_lbox[dim])
      r[dim] -= lbox[dim];
    else if (r[dim] < -half_lbox[dim])
      r[dim] += lbox[dim];
  }
}

inline void apply_pbc(double* r, int* b) { // periodic boundary conditions
  for (int dim = 0; dim < 3; ++dim) {
    if (r[dim] >= lbox[dim]) {
      r[dim] -= lbox[dim];
      b[dim]++;
    }
    else if (r[dim] < 0.0) {
      r[dim] += lbox[dim];
      b[dim]--;
    }
  }
}

void calc_lj_pair_force(double* r1, double* r2, double* f1, double* f2) {
  double r[3] = { r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2] };
	apply_min_img(r);
  const double rsq = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
  if (rsq < lj_rcut_sq) {
    const double s2 = lj_sigma * lj_sigma / rsq, s6 = s2 * s2 * s2, ff = -48.0 * lj_epsilon / rsq * s6 * (s6 - 0.5);
    for (int dim = 0; dim < 3; ++dim) {
      f1[dim] += ff * r[dim];
      f2[dim] -= ff * r[dim];
    }
  }
}

double lj_pair_energy(double* r1, double* r2) {
  double ener = 0.0, r[3] = { r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2] };
	apply_min_img(r);
  const double rsq = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
  if (rsq < lj_rcut_sq) {
    const double s2 = lj_sigma * lj_sigma / rsq, s6 = s2 * s2 * s2;
    ener = 4.0 * lj_epsilon * s6 * (s6 - 1.0) - lj_ucut;
  }
  return ener;
}

double lj_pair_vir(double* r1, double* r2) {
  double vir = 0.0, r[3] = { r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2] };
	apply_min_img(r);
  const double rsq = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
  if (rsq < lj_rcut_sq) {
    double s2 = lj_sigma * lj_sigma / rsq, s6 = s2 * s2 * s2;
    vir = 48.0 * lj_epsilon * s6 * (s6 - 0.5);
  }
  return vir;
}

double excess_pressure() { // excess pressure due to LJ interactions with tail correction
  double Pex = 0.0;
  for (int i = 0; i < n_particles; ++i) {
	  for (int k = 0; k < num_neigh[i]; ++k) {
		  int j = neigh_list[i][k];
		  Pex += lj_pair_vir(pos[i], pos[j]);
	  }
  }
  Pex = Pex / (3.0 * vol) + lj_ptail;
  return Pex;
}

void gen_neigh_list() { // generate Verlet list
  for (int i = 0; i < n_particles; ++i)
    num_neigh[i] = 0;
  for (int i = 0; i < n_particles; ++i)
    for (int j = i + 1; j < n_particles; ++j) {
      double r[3] = { pos[i][0] - pos[j][0], pos[i][1] - pos[j][1], pos[i][2] - pos[j][2] };
      apply_min_img(r);
      if (r[0] * r[0] + r[1] * r[1] + r[2] * r[2] <= nl_rcut_sq)
        neigh_list[i][num_neigh[i]++] = j;
    }
}

void update_neigh_list() {
  double max = 0.0, sec_max = 0.0, d2;
  for (int i = 0; i < n_particles; ++i) {
    d2 = dis[i][0] * dis[i][0] + dis[i][1] * dis[i][1] + dis[i][2] * dis[i][2];
    if (d2 > max) {
      sec_max = max;
      max = d2;
    }
    else if (d2 > sec_max)
      sec_max = d2;
  }
  if (sqrt(max) + sqrt(sec_max) > nl_rskin) {
    nl_update_count++;
    gen_neigh_list();
    for (int i = 0; i < n_particles; ++i)
      for (int dim = 0; dim < 3; ++dim)
        dis[i][dim] = 0.0;
  }
}

void calc_force() {
  for (int i = 0; i < n_particles; ++i) {
    for (int k = 0; k < num_neigh[i]; ++k) {
      int j = neigh_list[i][k];
      calc_lj_pair_force(pos[i], pos[j], fff[i], fff[j]);
    }
  }
}

void compute_energy_momentum(double& kin_ener, double& pot_ener, double* p) {
  kin_ener = 0.0;
  pot_ener = 0.0;
  for (int i = 0; i < n_particles; ++i) {
    for (int dim = 0; dim < 3; ++dim) {
      p[dim] += vel[i][dim];
      kin_ener += vel[i][dim] * vel[i][dim];
    }
    for (int k = 0; k < num_neigh[i]; ++k) {
      int j = neigh_list[i][k];
      pot_ener += lj_pair_energy(pos[i], pos[j]); 
    }
  }
  kin_ener *= 0.5;
  pot_ener += lj_utail * n_particles;
}

void integrate() {
  for (int i = 0; i < n_particles; ++i) {
		for (int dim = 0; dim < 3; ++dim) {
			vel[i][dim] += fff[i][dim] * md_delt_half; // velocity at t+dt/2
      fff[i][dim] = 0.0;
		  pos[i][dim] += vel[i][dim] * md_delt;
			dis[i][dim] += vel[i][dim] * md_delt;
		}
		apply_pbc(pos[i], pbc[i]);
  }
  update_neigh_list();
	// computing force field
	calc_force();
	for (int i = 0; i < n_particles; ++i)
		for (int dim = 0; dim < 3; ++dim)
			vel[i][dim] += fff[i][dim] * md_delt_half;  // velocity at t+dt
}

void corr_velocity() {
  double sum_v[3] = { 0.0, 0.0, 0.0 }, sum_v2 = 0.0;
  for (int i = 0; i < n_particles; ++i)
    for (int dim = 0; dim < 3; ++dim)
      sum_v[dim] += vel[i][dim];
  for (int dim = 0; dim < 3; ++dim)
    sum_v[dim] /= n_particles;
  for (int i = 0; i < n_particles; ++i)
    for (int dim = 0; dim < 3; ++dim) {
      vel[i][dim] -= sum_v[dim];
      sum_v2 += vel[i][dim] * vel[i][dim];
    }
  double kin_T = sum_v2 / (3.0 * (n_particles - 1) - 1), factor = sqrt(temp / kin_T);
  for (int i = 0; i < n_particles; ++i)
    for (int dim = 0; dim < 3; ++dim)
      vel[i][dim] *= factor;
}

void save_energy_momentum(int OPT) {
  static const int MAX = 1000000;
  static int K = 0;
  static double Ekin[MAX], Epot[MAX], p[MAX][3];
  switch (OPT) {
    case 1 : {
      if (K < MAX) {
        compute_energy_momentum(Ekin[K], Epot[K], p[K]);
        K++;
      }
    }
    break;
    case 2 : {
      std::ofstream out("ener_momentum.dat");
      out.precision(10);
      for (int k = 0; k < K; ++k)
        out << Ekin[k] << ' ' << Epot[k]  << ' ' << p[k][0] << ' ' << p[k][1] << ' ' << p[k][2] << '\n';
      out.close();
    }
    break;
  }
}

void save_pressure(int OPT) {
  static const int MAX = 1000000;
  static int K = 0;
  static double Pex[MAX];
  switch (OPT) {
    case 1 : {
      if (K < MAX)
        Pex[K++] = excess_pressure();
    }
    break;
    case 2 : {
      std::ofstream out("pressure.dat");
      out.precision(10);
      for (int k = 0; k < K; ++k)
        out << Pex[k] + rho * temp << '\n';
      out.close();
    }
    break;
  }
}

void calc_rdf(int OPT) {
  static const int MAX = 1000000;
  static int K = 0;
  static const int NBINS = 100;
  static const double rmin = 0.0, rmax = half_box, delt_r = (rmax - rmin) / NBINS;
  static double count[NBINS+1];
  switch (OPT) {
    case 1 : {
      if (K < MAX) {
        double r[3], dist;
        for (int i = 0; i < n_particles; ++i)
          for (int j = i + 1; j < n_particles; ++j) {
            for (int dim = 0; dim < 3; ++dim)
              r[dim] = pos[i][dim] - pos[j][dim];
            apply_min_img(r);
            dist = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
            int idx = int((dist - rmin) / delt_r + 0.0001);
            if (idx > NBINS) idx = NBINS;
            count[idx] += 2.0;
          }
        K++;
      }
    }
    break;
    case 2 : {
      std::ofstream out("rdf.dat");
      out.precision(10);
      for (int k = 0; k < NBINS; ++k) {
        double r = rmin + (k + 0.5) * delt_r;
        out << r << ' ' << count[k] / (4.0 * 3.14159265354 * r * r * delt_r * rho) / K / n_particles << '\n';
      }
      out.close();
    }
    break;
  }
}

void write_traj() {
  std::ofstream out("traj.xyz", std::ios::app);
  out.precision(10);
  out << n_particles << '\n';
  out << "LJ fluid\n";
  for (int n = 0; n < n_particles; ++n) {
    out << 'C' << "    " << pos[n][0] + pbc[n][0] * lbox[0]
               << "    " << pos[n][1] + pbc[n][1] * lbox[1] 
               << "    " << pos[n][2] + pbc[n][2] * lbox[2]
               << "    " << vel[n][0]
               << "    " << vel[n][1]
               << "    " << vel[n][2]
               << '\n';
  }
  out.close();
}

void save_syst_info() {
  std::ofstream out("syst.log");
  out << "Molecular dynamics simulations of Lennard-Jones fluid\n";
  out << "  Temperature: " << temp << '\n';
  out << "  " << n_particles << " particles in a cubic box of size " << box << " with the density " << rho << '\n';
  out << "  LJ potential parameters:\n";
  out << "    cutoff distance: " << lj_rcut << '\n';
  out << "    sigma: " << lj_sigma << '\n';
  out << "    epsilon: " << lj_epsilon << '\n';
  out << "  Run parameters: \n";
  out << "    relaxation steps: " << relax_steps << '\n';
  out << "    production steps: " << prod_steps << '\n';
  out << "    sampling steps: " << freq_sample << '\n';
  out << "    saving energy every " << freq_ener << " steps\n";
  out << "    scaling velocity every " << freq_scale << " steps\n";
  out << "    integration time step: " << md_delt << '\n';
  out << "  Verlet list parameters: \n";
  out << "    skin thickness: " << nl_rskin << '\n';
  out << "    cutoff radius: " << nl_rcut << '\n';
  out << "    maximum number of neighbors: " << max_num_neigh << '\n';
  out << "    total update times: " << nl_update_count << '\n';
  out.close();
}

void init_syst() {
  int seed = time(0);
  prng.init(seed);
  // read parameter from input file "parameter.txt"
  std::string sbuf;
  std::ifstream inputf("parameter.txt");
  inputf >> temp; getline(inputf, sbuf);
  inputf >> n_particles; getline(inputf, sbuf);
  inputf >> rho; getline(inputf, sbuf);
  inputf >> md_delt; getline(inputf, sbuf);
  inputf >> relax_steps; getline(inputf, sbuf);
  inputf >> prod_steps; getline(inputf, sbuf);
  inputf >> freq_sample; getline(inputf, sbuf);
  inputf >> freq_scale; getline(inputf, sbuf);
  inputf.close();
  // initialize other parameters and global variables
  beta = 1.0 / temp;
  pbc = new int*[n_particles];
  pos = new double*[n_particles];
  vel = new double*[n_particles];
  fff = new double*[n_particles];
  for (int n = 0; n < n_particles; ++n){
    pbc[n] = new int[3];
    pos[n] = new double[3];
    vel[n] = new double[3];
    fff[n] = new double[3];
    for (int dim = 0; dim < 3; ++dim)
      pbc[n][dim] = 0;
  }
  vol = n_particles / rho;
  box = pow(vol, 1.0 / 3.0);
  half_box = box * 0.5;
	lbox[0] = lbox[1] = lbox[2] = box;
  for (int dim = 0; dim < 3; ++dim)
    half_lbox[dim] = lbox[dim] * 0.5;
  lj_rcut = 2.5;
  lj_rcut_sq = lj_rcut * lj_rcut;
  lj_sigma = 1.0;
  lj_epsilon = 1.0;
  md_deltsq = md_delt * md_delt;
  md_delt_half = md_delt * 0.5;
  const double s1 = lj_sigma / lj_rcut, s2 = s1 * s1, s3 = s2 * s1, s6 = s3 * s3;
  lj_ucut = 4.0 * lj_epsilon * s6 * (s6 - 1.0);
  lj_utail = 8.0 / 3.0 * PI * rho * s3 * (s6 / 3.0 - 1.0);
  lj_ptail = 16.0 / 3.0 * PI * rho * rho * s3 * (s6 * 2.0 / 3.0 - 1.0);
  // put particles in a simple cubic lattice with density rho
  int n3 = int(pow((double)n_particles, 1.0 / 3.0)) + 1;
  double del = box / n3, dx = -del, dy, dz;
  int count = 0;
  for (int i = 0; i < n3; ++i) {
     dx += del;
     dy = -del;
     for (int j = 0; j < n3; ++j) {
        dy += del;
        dz = -del;
        for (int k = 0; k < n3; ++k) {
           dz += del;
           if (count < n_particles) {
              pos[count][0] = dx;
              pos[count][1] = dy;
              pos[count][2] = dz;
              count++;
           }
        }
     }
  }
  for (int i = 0; i < n_particles; ++i) 
    apply_pbc(pos[i], pbc[i]);

  // set velocities
  for (int i = 0; i < n_particles; ++i)
    for (int dim = 0; dim < 3; ++dim) {
      vel[i][dim] = prng.gen_open0_open1() - 0.5;
    }
  corr_velocity();
  // initialize Verlet list
  nl_rskin = 100.0 * sqrt(3.0 * temp / mass) * md_delt;
  nl_rcut = lj_rcut + nl_rskin;
  while (nl_rcut >= half_box) {
    nl_rskin *= 0.8;
    nl_rcut = lj_rcut + nl_rskin;
	}
	nl_rcut_sq = nl_rcut * nl_rcut;
  max_num_neigh = int((rho * PI * 2.0 / 3.0 * nl_rcut * nl_rcut * nl_rcut) * 4.0);
  num_neigh = new int[n_particles];
  neigh_list = new int*[n_particles]; // the zeroth index store the number of particles 
  for (int i = 0; i < n_particles; ++i)
    neigh_list[i] = new int[max_num_neigh];
  dis = new double*[n_particles];
  for (int i = 0; i < n_particles; ++i) {
    dis[i] = new double[3];
    for (int dim = 0; dim < 3; ++dim)
      dis[i][dim] = 0.0;
  }
  gen_neigh_list();
  for (int i = 0; i < n_particles; ++i)
    for (int dim = 0; dim < 3; ++dim)
      fff[i][dim] = 0.0;
  calc_force();
  save_energy_momentum(1);
}

int main() {
  init_syst();
  write_traj();
  for (int k = 1; k <= relax_steps; ++k) {
    integrate();
    if (k % freq_ener == 0) {
      save_energy_momentum(1);
    }
    if (k % freq_scale == 0) corr_velocity();  // scale velocities to get the desired temperature
    if (k % 1000 == 0) std::cout << '\r' << "relaxation run " << k << " / " << relax_steps << std::endl;
  }
  for (int k = 1; k <= prod_steps; ++k) {
    integrate();
    if (k % freq_ener == 0) {
      save_energy_momentum(1);
    }
    if (k % freq_sample == 0) {
      save_pressure(1);
      calc_rdf(1);
      write_traj();
    }
    if (k % 1000 == 0) std::cout << '\r' << "production run " << k << " / " << prod_steps << std::endl;
  }
  save_energy_momentum(2);
  save_pressure(2);
  calc_rdf(2);
  save_syst_info();
  return 0;
} 

