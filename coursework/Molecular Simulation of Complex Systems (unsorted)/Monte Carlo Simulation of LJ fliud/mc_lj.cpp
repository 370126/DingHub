/***************************************************
 NVT Monte Carlo simulations of Lennard-Jones fluid
 written by HU, Yi & HU, Jinglei 2016.03.09
***************************************************/
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include "prng.h"

#define PI 3.141592653589793238464

/* global variables in the simulations */
cPRNG<KISS, UNIFORM> prng; // random number generator;
int    n_particles; // total number of LJ particles
double temp, beta; // beta = 1.0 / kT
double rho, vol, box, half_box, lbox[3], half_lbox[3]; // density, volume and simulation box size
double lj_rcut, lj_rcut_sq, lj_sigma, lj_epsilon, lj_ucut, lj_utail, lj_ptail; // LJ parameters
double tot_ener; // total energy per particle
double **pos; // 2D array for particle positions (xi,yi,zi)
int    mc_freq_ener = 100, mc_relax_steps, mc_prod_steps, mc_freq_sample; // mc relaxation and production steps
double mc_max_dis = 1.0, mc_ntry = 0.0, mc_nsucc = 0.0;

/* functions */
void apply_min_img(double* r) { // minimum image convention
  for (int dim = 0; dim < 3; ++dim) {
    if (r[dim] > half_lbox[dim])
      r[dim] -= lbox[dim];
    else if (r[dim] < -half_lbox[dim])
      r[dim] += lbox[dim];
  }
}

void apply_pbc(double* r) { // periodic boundary conditions
  for (int dim = 0; dim < 3; ++dim) {
    if (r[dim] >= lbox[dim])
      r[dim] -= lbox[dim];
    else if (r[dim] < 0.0)
      r[dim] += lbox[dim];
  }
} // 0--Lx, 0--Ly, 0--Lz --> -Lx/2 -- +Lx/2, -Ly/2 -- +Ly/2, -Lz/2 -- +Lz/2

double lj_pair_ener(double* r1, double* r2) {
  double r[3] = { r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2] };
  apply_min_img(r);
  double ener = 0.0, rsq = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
  if (rsq < lj_rcut_sq) {
    double s2 = lj_sigma * lj_sigma / rsq, s6 = s2 * s2 * s2;
    ener = 4.0 * lj_epsilon * s6 * (s6 - 1.0) - lj_ucut;
  }
  return ener;
}

double lj_pair_vir(double* r1, double* r2) {
  double r[3] = { r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2] };
  apply_min_img(r);
  double vir = 0.0, rsq = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
  if (rsq < lj_rcut_sq) {
    double s2 = lj_sigma * lj_sigma / rsq, s6 = s2 * s2 * s2;
    vir = 48.0 * lj_epsilon * s6 * (s6 - 0.5);
  }
  return vir;
}

double energy_per_particle()
 { // total energy per particle with tail correction
  double Etot = 0.0;
  for (int i = 0; i < n_particles - 1; ++i) {
    for (int j = i + 1; j < n_particles; ++j) {
      Etot += lj_pair_ener(pos[i], pos[j]);
    }
  }
  Etot = Etot / n_particles + lj_utail + lj_ucut;
  return Etot;
}

double excess_pressure() { // excess pressure due to LJ interactions with tail correction
  double Pex = 0.0;
  for (int i = 0; i < n_particles - 1; ++i)
    for (int j = i + 1; j < n_particles; ++j)
      Pex += lj_pair_vir(pos[i], pos[j]);
  Pex = Pex / (3.0 * vol) + lj_ptail;
  return Pex;
}

void mc_sweep() { // trying to move all particles once during one MC sweep or step
  for (int n = 0; n < n_particles; ++n) {
    const int i = int(n_particles * prng.gen_open0_open1()); // randomly choose a particle i
    double ener_old = 0.0, ener_new = 0.0, ri_new[3] = { pos[i][0] + (prng.gen_open0_open1() - 0.5) * mc_max_dis, \
                                                         pos[i][1] + (prng.gen_open0_open1() - 0.5) * mc_max_dis, \
                                                         pos[i][2] + (prng.gen_open0_open1() - 0.5) * mc_max_dis };
    //double ener_old = 0.0, ener_new = 0.0, ri_new[3] = { pos[i][0] + (prng.gen_open0_open1() - 0.) * mc_max_dis, \
                                                         pos[i][1] + (prng.gen_open0_open1() - 0.) * mc_max_dis, \
                                                         pos[i][2] + (prng.gen_open0_open1() - 0.) * mc_max_dis };
    apply_pbc(ri_new);
    for (int j = 0; j < n_particles; ++j) {
      if (j != i) {
        ener_old += lj_pair_ener(pos[i], pos[j]);
        ener_new += lj_pair_ener(ri_new, pos[j]);
      }
    }
    double dE = ener_new - ener_old;
    mc_ntry++;
    if (dE <= 0.0 || exp(-beta * dE) > prng.gen_close0_open1()) { // accept the trial move
      for (int dim = 0; dim < 3; ++dim)
        pos[i][dim] = ri_new[dim];
      tot_ener += dE / n_particles; //! tot_ener: total energy per particle
      mc_nsucc++;
    } 
  }
}

void adjust_mc_move(double acc_ratio) { // determine max. displacement of trial move such that the acc ratio is between 30% and 50%
  if (acc_ratio < 0.3) {
    mc_max_dis *= 0.9;
  } else if (acc_ratio > 0.5) {
    mc_max_dis *= 1.2;
    if (mc_max_dis >= half_box)
      mc_max_dis = half_box;
  }
}

void save_energy(int OPT) {
  static const int MAX = 1000000;
  static int K = 0;
  static double ener_1[MAX], ener_2[MAX];
  switch (OPT) {
    case 1 : {
      if (K < MAX) {
        ener_1[K] = tot_ener;
        ener_2[K] = energy_per_particle();
        K++;
      }
    }
    break;
    case 2 : {
      std::ofstream out("tot_ener.dat");
      out.precision(10);
      for (int k = 0; k < K; ++k)
        out << ener_1[k] << ' ' << ener_2[k] << '\n';
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

void write_traj() {
  std::ofstream out("traj.xyz", std::ios::app);
  out.precision(10);
  out << n_particles << '\n';
  out << "LJ fluid\n";
  for (int n = 0; n < n_particles; ++n) {
    out << 'C'
        << "    " << pos[n][0]
        << "    " << pos[n][1]
        << "    " << pos[n][2]
        << '\n';
  }
  out.close();
}

void init_syst() {
  // read parameter from input file "parameter.txt"
  std::string sbuf;
  std::ifstream inputf("parameter.txt");
  inputf >> temp; getline(inputf, sbuf);
  inputf >> n_particles; getline(inputf, sbuf);
  inputf >> rho; getline(inputf, sbuf);
  inputf >> lj_rcut; getline(inputf, sbuf);
  inputf >> mc_relax_steps; getline(inputf, sbuf);
  inputf >> mc_prod_steps; getline(inputf, sbuf);
  inputf >> mc_freq_sample; getline(inputf, sbuf);
  inputf.close();
  // initialize other parameters and global variables
  beta = 1.0 / temp;
  pos = new double*[n_particles];
  for (int n = 0; n < n_particles; ++n)
    pos[n] = new double[3];
  vol = n_particles / rho;
  box = pow(vol, 1.0 / 3.0);
  half_box = box * 0.5;
	lbox[0] = lbox[1] = lbox[2] = box;
  for (int dim = 0; dim < 3; ++dim)
    half_lbox[dim] = lbox[dim] * 0.5;
  if (lj_rcut > half_box) lj_rcut = half_box;
  lj_rcut_sq = lj_rcut * lj_rcut;
  lj_sigma = 1.0;
  lj_epsilon = 1.0;
  double s1 = lj_sigma / lj_rcut;
  double s2 = s1*s1;
  double s3 = s2*s1;
  double s6 = s3*s3;
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
  // fold positions back to the central box
  for (int n = 0; n < n_particles; ++n)
    apply_pbc(pos[n]);
  tot_ener = energy_per_particle();
}

void save_syst_info() {
  std::ofstream out("sys_info.log");
  out << "Monte Carlo simulations of Lennard-Jones fluid\n";
  out << "  Temperature: " << temp << '\n';
  out << "  " << n_particles << " particles in a cubic box of size " << box << " with the density " << rho << '\n';
  out << "  LJ potential parameters:\n";
  out << "    cutoff distance: " << lj_rcut << '\n';
  out << "    sigma: " << lj_sigma << '\n';
  out << "    epsilon: " << lj_epsilon << '\n';
  out << "  MC parameters: \n";
  out << "    relaxation steps: " << mc_relax_steps << '\n';
  out << "    production steps: " << mc_prod_steps << '\n';
  out << "    sampling steps: " << mc_freq_sample << '\n';
  out << "    maximum displacement: " << mc_max_dis << '\n';
  out << "    acceptance ratio (%): " << mc_nsucc / mc_ntry * 100.0 << '\n';
  out.close();
}

int main() {
  int seed = time(0);
  prng.init(seed);
  init_syst();
  for (int k = 0; k < mc_relax_steps / 10; ++k) {
    mc_sweep();
    if (k % 5 == 0) {
      double acc_ratio = mc_nsucc / mc_ntry;
      adjust_mc_move(acc_ratio);
      mc_ntry = 0.0;
      mc_nsucc = 0.0;
    }
  }
  //save_energy(1);
  for (int k = 0; k < mc_relax_steps; ++k) { // relaxation run towards equilibrium
    //if (k % mc_freq_ener == 0) save_energy(1);
    mc_sweep();
  }
  for (int k = 0; k < mc_prod_steps; ++k) { // production run
    mc_sweep();
    if (k % mc_freq_ener == 0) save_energy(1);
    if (k % mc_freq_sample == 0) {
      save_pressure(1);
      write_traj();
    }
  }
  save_energy(2);
  save_pressure(2);
  save_syst_info();
  return 0;
} 

