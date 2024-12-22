/****************************************************************************************
  Dissipative Particle Dynamics simulations of lipids (DMPC) bilayer
  Copyright (c) 2017.06 Jinglei Hu < hujinglei _at_ nju.edu.cn >
  Macro usage:
    _OPENMP: openmp parallelization
  NOTE: Refer to Andrea Grafmueller's PhD thesis (P41, Eq. 3.17)
        for the calculation of stress tensor!
****************************************************************************************/
/* header files */
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory.h>
#ifdef _OPENMP
#include <omp.h>                  // OpenMP header file
#endif
#include <sstream>
#include <string>
#include "prng.h"                 // Random number generator library
#include "timer.h"                // Timer Class

#ifdef _OPENMP
#define N_THREADS                 6
#else
#define N_THREADS                 1
#endif
#define FREQ_SORT                 25
#define RESTART                   false
#define RELAX                     true

/* Constants */
#define HALF_SQRT3                0.8660254037844386

/* Parameters */
#define k_BT                      1.0
#define RHO                       3.0
#define BOXLX                     20.0
#define BOXLY                     20.0
#define BOXLZ                     20.0
#define VOL                       (BOXLX * BOXLY * BOXLZ)
#define INV_BOXLX                 (1.0 / BOXLX)
#define INV_BOXLY                 (1.0 / BOXLY)
#define INV_BOXLZ                 (1.0 / BOXLZ)
#define N_BEADS                   (int)(VOL * RHO + 0.5)
#define DPD_R0                    1.0
#define DPD_R0_SQ                 (DPD_R0 * DPD_R0)
#define CELL_RC                   DPD_R0
#define N_CELLS_X                 (int)(BOXLX / CELL_RC + 1.0e-09)
#define N_CELLS_Y                 (int)(BOXLY / CELL_RC + 1.0e-09)
#define N_CELLS_Z                 (int)(BOXLZ / CELL_RC + 1.0e-09)
#define N_CELLS                   (N_CELLS_X * N_CELLS_Y * N_CELLS_Z)
#define INV_CELL_LX               (N_CELLS_X * INV_BOXLX)
#define INV_CELL_LY               (N_CELLS_Y * INV_BOXLY)
#define INV_CELL_LZ               (N_CELLS_Z * INV_BOXLZ)
#define DPD_DT                    0.03
#define HALF_DPD_DT               (DPD_DT * 0.5)
#define INIT_CYCLE                10000
#define PROD_CYCLE                50000
#define FREQ_RESCALE              50
#define FREQ_SAMPLE               10
#define N_CPT                     10
#define FREQ_CPT                  (PROD_CYCLE / N_CPT)
#define N_SAMPLES                 (PROD_CYCLE / FREQ_SAMPLE)
#define FREQ_LOG                  1000

/* Lipid */
#define AREA_PER_LIPID            1.2
#define N_HEADS_PER_LIPID         3
#define N_BEADS_PER_TAIL          4
#define N_BEADS_PER_LIPID         (N_HEADS_PER_LIPID + N_BEADS_PER_TAIL * 2)
#define K_BOND_LIPID              128.0   // bond spring constant
#define L0_BOND_LIPID             0.5     // and equilibrium length
#define K_ANGLE_LIPID             15.0    // angle potential strength
#define PHI0_CCC                  0.0     // and preferred angle
#define N_LIPIDS                  (2 * int(BOXLX * BOXLY / AREA_PER_LIPID + 1.0e-09))
#define N_LIPID_BEADS             (N_LIPIDS * N_BEADS_PER_LIPID)
#define N_SOLVENT                 (N_BEADS - N_LIPID_BEADS)

/* Global variables */
cPRNG<MT19937, NORMAL::ZIGGURAT> prng[N_THREADS];
const double BOXL[3] = { BOXLX, BOXLY, BOXLZ }, HALF_BOXL[3] = { BOXLX * 0.5, BOXLY * 0.5, BOXLZ * 0.5 };
double       dpd_f_par[3][3][3] = { 0.0 }; // dpd force parameters a_ij, gamma_ij & sigma_ij = sqrt(2 k_BT gamma_ij/ Dt)
int          bonded_pair[N_LIPIDS*10][2], bonded_triple[N_LIPIDS*6][3];
bool         at_boundary[N_CELLS][13]; // true if the neighboring cell is at boundary, or false otherwise
int          cell_list[N_CELLS][13], head_of_cell[N_CELLS], cell_idx[N_BEADS], linked_list[N_BEADS], idx_of_occur[N_SOLVENT];
int          box[N_LIPID_BEADS][3] = { 0 };
double       pos[N_BEADS][3], vel[N_BEADS][3], fff[N_THREADS][N_BEADS+6][3] = { 0.0 }, _pos[N_SOLVENT][3], _vel[N_SOLVENT][3], _fff[N_SOLVENT][3];
double       t_force = 0.0, t_list = 0.0, t_sort = 0.0, t_integrate = 0.0, Ekin = 0.0, Epot[N_THREADS][18] = { 0.0 };
#ifdef TENSION
double       stress_tensor[N_THREADS][9+16] = { 0.0 };
#endif

/* Functions */
void apply_min_img(double* r) {
  for (int dim = 0; dim < 3; ++dim) {
    if (r[dim] > HALF_BOXL[dim])
      r[dim] -= BOXL[dim];
    else if (r[dim] < -HALF_BOXL[dim])
      r[dim] += BOXL[dim];
  }
}

void apply_pbc(double* r) {
  for (int dim = 0; dim < 3; ++dim) {
    if (r[dim] >= BOXL[dim])
      r[dim] -= BOXL[dim];
    else if (r[dim] < 0.0)
      r[dim] += BOXL[dim];
  }
}

void apply_pbc(double* r, int* b) {
  for (int dim = 0; dim < 3; ++dim) {
    if (r[dim] >= BOXL[dim]) {
      r[dim] -= BOXL[dim];
      b[dim] ++;
    }
    else if (r[dim] < 0.0) {
      r[dim] += BOXL[dim];
      b[dim] --;
    }
  }
}

void save_check_point() {
  std::ofstream out("dpd.cpt");
  out.precision(6);
  for (int n = 0; n < N_LIPID_BEADS; ++n) // solute
    out << std::setw(10) << pos[n][0] << ' ' 
        << std::setw(10) << pos[n][1] << ' '
        << std::setw(10) << pos[n][2] << ' '
        << std::setw(10) << vel[n][0] << ' '
        << std::setw(10) << vel[n][1] << ' '
        << std::setw(10) << vel[n][2] << ' '
        << std::setw(8) << box[n][0] << ' '
        << std::setw(8) << box[n][1] << ' '
        << std::setw(8) << box[n][2] << '\n';
  for (int n = N_LIPID_BEADS; n < N_BEADS; ++n) // solvent
    out << std::setw(10) << pos[n][0] << ' ' 
        << std::setw(10) << pos[n][1] << ' '
        << std::setw(10) << pos[n][2] << ' '
        << std::setw(10) << vel[n][0] << ' '
        << std::setw(10) << vel[n][1] << ' '
        << std::setw(10) << vel[n][2] << '\n';
  out.close(); 
}

void write_log() {
  std::ofstream out("sys.log");
  out << "---===System===---\n";
  out << "Bilayer membrane in water\n";
  out << "Temperature: " << k_BT << '\n';
  out << "Box: " << BOXLX << " * " << BOXLY << " * " << BOXLZ << '\n';
  out << "RHO: " << RHO << '\n';
  out << "Total Number of Beads: " << N_BEADS << '\n';
  out << " Water beads: " << N_SOLVENT << '\n';
  out << " Lipid beads: " << N_LIPID_BEADS << '\n';
  out << "Number of Lipid: " << N_LIPIDS << '\n';
  out << "Number of Beads per Lipid: " << N_BEADS_PER_LIPID << '\n';
  out << "Concentration of Lipid: " << N_LIPID_BEADS / VOL << '\n';
  out << "Predefined projected area per Lipid: " << AREA_PER_LIPID << '\n';
  out << "Real projected area per Lipid: " << 2.0 * BOXLX * BOXLY / N_LIPIDS << '\n';
  out << "Lipid model parameters:\n";
  out << " spring constant: " << K_BOND_LIPID << '\n';
  out << " preferred bond length: " << L0_BOND_LIPID << '\n';
  out << " angle potential strength: " << K_ANGLE_LIPID << '\n';
  out << " preferred bond angle: " << PHI0_CCC << '\n';
  out << "\n---===DPD===---\n";
  #ifdef _OPENMP
  out << "Parallel run with " << N_THREADS << " CPU cores\n";
  #else
  out << "Serial run\n";
  #endif
  out << "Sort solvent particles every " << FREQ_SORT << " time steps\n";
  out << "Start from checkpoint: " << (RESTART ? "Yes" : "No") << '\n';
  out << "Time step: " << DPD_DT << '\n';
  out << "Number of cycles:\n";
  if (RELAX) out << " Relaxation run: " << INIT_CYCLE << '\n';
  out << " Production run: " << PROD_CYCLE << '\n';
  out << "Sampling frequency: " << FREQ_SAMPLE << '\n';
  out << "Beads type: ";
  const std::string sBeadsType[] = { "W", "H", "C" };
  for (int i = 0; i < 3; ++i)
    out << i << " - " << sBeadsType[i] << ' ';
  out << "\nDPD parameters:\n";
  out << " a_ij(g_ij)\t";
  for (int i = 0; i < 3; ++i)
    out << sBeadsType[i] << '\t';
  out << '\n';
  for (int i = 0; i < 3; ++i) {
    out << ' ' << sBeadsType[i] << "\t\t";
    for (int j = 0; j < 3; ++j)
      out << dpd_f_par[i][j][0] << '(' << dpd_f_par[i][j][1] << ")\t";
    out << '\n';
  }
  out << "Running...\n";
  out.close();
}

void write_traj(int OPT = 1) {
  switch (OPT) {
    case 0: {
        std::ofstream out("mol.psf");
        out << "PSF\n";
        out << std::setw(8) << std::right << N_LIPID_BEADS << " !NATOM\n";
        for (int i = 0; i < N_LIPID_BEADS; ++i)
          out << std::setw(8) << std::right << i+1
            << std::setw(2) << std::right << 'X'
            << std::setw(5) << std::right << 0
            << std::setw(7) << std::right << (i % N_BEADS_PER_LIPID <= 2 ? 'O' : 'C')
            << std::setw(3) << std::right << (i % N_BEADS_PER_LIPID <= 2 ? 'O' : 'C')
            << std::setw(3) << std::right << (i % N_BEADS_PER_LIPID <= 2 ? 'O' : 'C')
            << std::setw(21) << ' ' << '\n';
        out << std::setw(8) << std::right << N_LIPIDS * 10 << " !NBOND\n";
        int kk = 0;
        for (int k = 0; k < (N_LIPIDS * 10) / 4; ++k) {
          out << std::setw(8) << std::right << bonded_pair[kk][0]  + 1 << std::setw(8) << std::right << bonded_pair[kk][1] + 1
              << std::setw(8) << std::right << bonded_pair[kk+1][0]  + 1 << std::setw(8) << std::right << bonded_pair[kk+1][1] + 1
              << std::setw(8) << std::right << bonded_pair[kk+2][0]  + 1 << std::setw(8) << std::right << bonded_pair[kk+2][1] + 1
              << std::setw(8) << std::right << bonded_pair[kk+3][0]  + 1 << std::setw(8) << std::right << bonded_pair[kk+3][1] + 1
              << '\n';
          kk += 4;
        }
        for (int k = ((N_LIPIDS * 10) / 4) * 4; k < N_LIPIDS * 10; ++k) {
          out << std::setw(8) << std::right << bonded_pair[kk][0]  + 1 << std::setw(8) << std::right << bonded_pair[kk][1] + 1;
          kk++;
        }
        out.close();
      }
      break;
    case 1: {
        std::ofstream out("traj.xyz", std::ios::app);
        out.precision(10);
        out << N_LIPID_BEADS << '\n';
        out << "dpd lipid bilayer\n";
        for (int n = 0; n < N_LIPID_BEADS; ++n) {
          out << 'X'
              << "    " << pos[n][0] << "    " << pos[n][1] << "    " << pos[n][2]
              << "    " << box[n][0] << "    " << box[n][1] << "    " << box[n][2]
              << '\n';
        }
        out.close();
      }
      break;
  }
}

void gen_cell_list() {
  const int dcell[13][3] = { {1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {1, 1, 0}, 
                             {-1, -1, 1}, {0, -1, 1}, {1, -1, 1}, {-1, 0, 1}, {0, 0, 1}, {1, 0, 1}, {-1, 1, 1}, {0, 1, 1}, {1, 1, 1} };
  for (int ix = 0; ix < N_CELLS_X; ++ix)
    for (int iy = 0; iy < N_CELLS_Y; ++iy)
      for (int iz = 0; iz < N_CELLS_Z; ++iz) {
        const int celli = ix + (iy + iz * N_CELLS_Y) * N_CELLS_X;
        for (int n = 0; n < 13; ++n) {
          int jx = ix + dcell[n][0], jy = iy + dcell[n][1], jz = iz + dcell[n][2];
          at_boundary[celli][n] = false;
          if (jx < 0) { jx += N_CELLS_X; at_boundary[celli][n] = true; }
            else if (jx >= N_CELLS_X) { jx -= N_CELLS_X; at_boundary[celli][n] = true; }
          if (jy < 0) { jy += N_CELLS_Y; at_boundary[celli][n] = true; }
            else if (jy >= N_CELLS_Y) { jy -= N_CELLS_Y; at_boundary[celli][n] = true; }
          if (jz < 0) { jz += N_CELLS_Z; at_boundary[celli][n] = true; }
            else if (jz >= N_CELLS_Z) { jz -= N_CELLS_Z; at_boundary[celli][n] = true; }
          cell_list[celli][n] = jx + (jy + jz * N_CELLS_Y) * N_CELLS_X;
        }
      }
} // the indices of 13 neighbor cells for each cell

void gen_linked_list() {
  cTimer timer;
  for (int cell = 0; cell < N_CELLS; ++cell)
    head_of_cell[cell] = -1;
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static, N_BEADS / N_THREADS) num_threads(N_THREADS)
  #endif
  for (int n = 0; n < N_BEADS; ++n)
    cell_idx[n] = int(pos[n][0] * INV_CELL_LX) + (int(pos[n][1] * INV_CELL_LY) + int(pos[n][2] * INV_CELL_LZ) * N_CELLS_Y) * N_CELLS_X;
  for (int n = 0; n < N_BEADS; ++n) {
    const int cell = cell_idx[n];
    linked_list[n] = head_of_cell[cell];
    head_of_cell[cell] = n;
  }
  t_list += timer.elapsed();
} // generate cell-linked list

void sort_solvent() {
  cTimer timer;
  for (int cell = 0; cell < N_CELLS; ++cell)
    head_of_cell[cell] = -1;
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static, N_SOLVENT / N_THREADS) num_threads(N_THREADS)
  #endif
  for (int n = N_LIPID_BEADS; n < N_BEADS; ++n)
    cell_idx[n] = int(pos[n][0] * INV_CELL_LX) + (int(pos[n][1] * INV_CELL_LY) + int(pos[n][2] * INV_CELL_LZ) * N_CELLS_Y) * N_CELLS_X;
  for (int n = N_LIPID_BEADS; n < N_BEADS; ++n) {
    const int cell = cell_idx[n];
    linked_list[n] = head_of_cell[cell];
    head_of_cell[cell] = n;
  }
  int counts = 0;
  for (int cell = 0; cell < N_CELLS; ++cell) {
    int n = head_of_cell[cell];
    while (n > -1) {
      idx_of_occur[counts++] = n;
      n = linked_list[n];
    }
  }
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static, N_SOLVENT / N_THREADS)  num_threads(N_THREADS)
  #endif
  for (int n = 0; n < counts; ++n) { // counts = N_SOLVENT
    const int k = idx_of_occur[n];
    for (int dim = 0; dim < 3; ++dim) {
      _pos[n][dim] = pos[k][dim];
      _vel[n][dim] = vel[k][dim];
      _fff[n][dim] = fff[0][k][dim];
    }
  }
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static, N_SOLVENT / N_THREADS)  num_threads(N_THREADS)
  #endif
  for (int n = 0; n < N_SOLVENT; ++n) {
    const int offset = N_LIPID_BEADS;
    for (int dim = 0; dim < 3; ++dim) {
      pos[n+offset][dim] = _pos[n][dim];
      vel[n+offset][dim] = _vel[n][dim];
      fff[0][n+offset][dim] = _fff[n][dim];
    }
  }
  t_sort += timer.elapsed();
} // to ensure solvent particles close in simulation box are stored close in memory

int bead_type(int n) {
  if (n >= N_LIPID_BEADS)
    return 0; // water
  else if (n % N_BEADS_PER_LIPID > 2)
    return 2; // tail
  else
    return 1; // head
} //! NOTE: in pos, vel and fff arrays, first lipid beads, then water beads

void calc_dpd_force(int i, int j, int TID = 0, bool MinImg = true) {
  double rij[] = { pos[i][0] - pos[j][0], pos[i][1] - pos[j][1], pos[i][2] - pos[j][2] };
  if (MinImg) apply_min_img(rij);
  const double rsq = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];
  if (rsq < DPD_R0_SQ) {
    const double r = sqrt(rsq), w = 1.0 / r - 1.0;
    const int ti = bead_type(i), tj = bead_type(j);
    const double a_ij = dpd_f_par[ti][tj][0], g_ij = dpd_f_par[ti][tj][1], s_ij = dpd_f_par[ti][tj][2];
    Epot[TID][0] += 0.5 * a_ij * (1.0 + rsq - r - r);
    const double vij[] = { vel[i][0] - vel[j][0], vel[i][1] - vel[j][1], vel[i][2] - vel[j][2] }, \
                 rv = rij[0] * vij[0] + rij[1] * vij[1] + rij[2] * vij[2], \
                 coeff = (a_ij - g_ij * w * rv + s_ij * prng[TID].gen_gauss()) * w, \
                 ff[] = { rij[0] * coeff, rij[1] * coeff, rij[2] * coeff };
    for (int dim = 0; dim < 3; ++dim) {
      fff[TID][i][dim] += ff[dim];
      fff[TID][j][dim] -= ff[dim];
    }
    #ifdef TENSION
    for (int dim1 = 0; dim1 < 3; ++dim1)
      for (int dim2 = 0; dim2 < 3; ++dim2)
        stress_tensor[TID][3*dim1+dim2] -= ff[dim1] * rij[dim2]; // It's also fine to include conservative force (a_ij * w) only.
    #endif
  }
}

void calc_bond_force(int i, int j, double Kb, double L0, int TID = 0) {
  double rij[] = { pos[i][0] - pos[j][0], pos[i][1] - pos[j][1], pos[i][2] - pos[j][2] };
  apply_min_img(rij);
  const double r = sqrt(rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]), d = r - L0, coeff = Kb * (L0 / r - 1.0);
  Epot[TID][0] += 0.5 * Kb * d * d;
  const double ff[] = { rij[0] * coeff, rij[1] * coeff, rij[2] * coeff };
  for (int dim = 0; dim < 3; ++dim) {
    fff[TID][i][dim] += ff[dim];
    fff[TID][j][dim] -= ff[dim];
  }
  #ifdef TENSION
  for (int dim1 = 0; dim1 < 3; ++dim1)
    for (int dim2 = 0; dim2 < 3; ++dim2)
      stress_tensor[TID][3*dim1+dim2] -= ff[dim1] * rij[dim2];
  #endif
}

void calc_angle_force(int i, int j, int k, double k_theta, int TID = 0) {
  double rji[] = { pos[j][0] - pos[i][0], pos[j][1] - pos[i][1], pos[j][2] - pos[i][2] }, \
         rkj[] = { pos[k][0] - pos[j][0], pos[k][1] - pos[j][1], pos[k][2] - pos[j][2] };
  apply_min_img(rji);
  apply_min_img(rkj);
  const double inv_rjisq = 1.0 / (rji[0] * rji[0] + rji[1] * rji[1] + rji[2] * rji[2]), \
               inv_rkjsq = 1.0 / (rkj[0] * rkj[0] + rkj[1] * rkj[1] + rkj[2] * rkj[2]), \
               inv_rji_rkj = sqrt(inv_rjisq * inv_rkjsq), \
               rji_dot_rkj = rji[0] * rkj[0] + rji[1] * rkj[1] + rji[2] * rkj[2];
  Epot[TID][0] += k_theta * (1.0 - inv_rji_rkj * rji_dot_rkj);
  const double coeff = k_theta * inv_rji_rkj;
  double fi[3], fk[3];
  for (int dim = 0; dim < 3; ++dim) {
    fi[dim] = coeff * (rji_dot_rkj * inv_rjisq * rji[dim] - rkj[dim]);
    fk[dim] = coeff * (rji[dim] - rji_dot_rkj * inv_rkjsq * rkj[dim]);
  }
  double fj[] = { -fi[0]-fk[0], -fi[1]-fk[1], -fi[2]-fk[2] };
  for (int dim = 0; dim < 3; ++dim) {
    fff[TID][i][dim] += fi[dim];
    fff[TID][j][dim] += fj[dim];
    fff[TID][k][dim] += fk[dim];
  }
  #ifdef TENSION
  double rik[] = { pos[i][0]-pos[k][0], pos[i][1]-pos[k][1], pos[i][2]-pos[k][2] };
  apply_min_img(rik);
  for (int dim1 = 0; dim1 < 3; ++dim1)
    for (int dim2 = 0; dim2 < 3; ++dim2)
      stress_tensor[TID][3*dim1+dim2] -= 0.3333333333333333 * ((fj[dim1]-fi[dim1])*rji[dim2] + (fk[dim1]-fj[dim1])*rkj[dim2] + (fi[dim1]-fk[dim1])*rik[dim2]);
  #endif
}

void compute_forces() {
  gen_linked_list();
  cTimer timer;
  for (int t = 0; t < N_THREADS; ++t) {
    Epot[t][0] = 0.0;
    #ifdef TENSION
    for (int i = 0; i < 9; ++i)
      stress_tensor[t][i] = 0.0;
    #endif
    for (int k = 0; k < N_BEADS; ++k)
      for (int dim = 0; dim < 3; ++dim)
        fff[t][k][dim] = 0.0;
  }
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static, N_LIPID_BEADS * 10 /N_THREADS) num_threads(N_THREADS)
  #endif
  for (int n = 0; n < N_LIPIDS * 10; n++) {
    #ifdef _OPENMP
    int TID = omp_get_thread_num();
    #else
    int TID = 0;
    #endif
    calc_bond_force(bonded_pair[n][0], bonded_pair[n][1], K_BOND_LIPID, L0_BOND_LIPID, TID);
  } // bond force
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static, N_LIPID_BEADS * 6 /N_THREADS) num_threads(N_THREADS)
  #endif
  for (int n = 0; n < N_LIPIDS * 6; n++) {
    #ifdef _OPENMP
    int TID = omp_get_thread_num();
    #else
    int TID = 0;
    #endif
    calc_angle_force(bonded_triple[n][0], bonded_triple[n][1], bonded_triple[n][2], K_ANGLE_LIPID, TID);
  } // angle force
  #ifdef _OPENMP
  #pragma omp parallel for schedule(guided) num_threads(N_THREADS)
  #endif
  for (int celli = 0; celli < N_CELLS; ++celli) {
    #ifdef _OPENMP
    int TID = omp_get_thread_num();
    #else
    int TID = 0;
    #endif
    int i = head_of_cell[celli];
    while (i > -1) {
      int j = linked_list[i];
      while (j > -1) {
        calc_dpd_force(i, j, TID, false);
        j = linked_list[j];
      } // loop over beads in the same cell
      for (int k = 0; k < 13; ++k) {
        bool flag = at_boundary[celli][k];
        int cellj = cell_list[celli][k];
        j = head_of_cell[cellj];
        while (j > -1) {
          calc_dpd_force(i, j, TID, flag);
          j = linked_list[j];
        }
      } // loop over beads in the 13 neighboring cells
      i = linked_list[i];
    }
  } // nonbonded force
  // sum up
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static, N_BEADS/N_THREADS) num_threads(N_THREADS)
  for (int n = 0; n < N_BEADS; ++n)
    for (int t = 1; t < N_THREADS; ++t)
      for (int dim = 0; dim < 3; ++dim)
        fff[0][n][dim] += fff[t][n][dim];
  for (int t = 1; t < N_THREADS; ++t)
    Epot[0][0] += Epot[t][0];
  #endif
  t_force += timer.elapsed();
}

void integrate() {
  cTimer timer;
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static, N_BEADS/N_THREADS) num_threads(N_THREADS)
  #endif
  for (int n = 0; n < N_BEADS; ++n) {
    for (int dim = 0; dim < 3; ++dim) {
      vel[n][dim] += fff[0][n][dim] * HALF_DPD_DT;  // new f dependes on v, not precise but very cheap
      pos[n][dim] += vel[n][dim] * DPD_DT;
    }
    if (n < N_LIPID_BEADS)
      apply_pbc(pos[n], box[n]);
    else
      apply_pbc(pos[n]);
  }
  t_integrate += timer.elapsed();
  compute_forces();
  timer.restart();
  Ekin = 0.0;
  #ifdef _OPENMP
  #pragma omp parallel for reduction(+:Ekin) schedule(static, N_BEADS/N_THREADS) num_threads(N_THREADS)
  #endif
  for (int n = 0; n < N_BEADS; ++n) {
    for (int dim = 0; dim < 3; ++dim)
      vel[n][dim] += fff[0][n][dim] * HALF_DPD_DT;
    Ekin += vel[n][0] * vel[n][0] + vel[n][1] * vel[n][1] + vel[n][2] * vel[n][2];
  }
  t_integrate += timer.elapsed();
} // Velocity-Verlet algorithm ($\lambda = 0.5$)

void rescale_Tkin() {
  double coeff = sqrt(3.0 * (N_BEADS - 1) * k_BT / Ekin);
  for (int i = 0; i < N_BEADS; ++i)
    for (int dim = 0; dim < 3; ++dim)
      vel[i][dim] *= coeff;
  //! Zero total momentum
  double Impulse[] = { 0.0, 0.0, 0.0 };
  for (int i = 0; i < N_BEADS; ++i)
    for (int dim = 0; dim < 3; ++dim)
      Impulse[dim] += vel[i][dim];
  for (int i = 0; i < N_BEADS; ++i)
    for (int dim = 0; dim < 3; ++dim)
      vel[i][dim] -= Impulse[dim] / (double)N_BEADS;
}

void gen_init_conf() {
  const double thickness = (N_BEADS_PER_TAIL + 1) * L0_BOND_LIPID, z[] = { BOXLZ * 0.5 - thickness, BOXLZ * 0.5 + thickness };
  const int num_lipid[] = {  N_LIPIDS / 2, N_LIPIDS -  N_LIPIDS / 2 };
  const int ngridx = (int)(sqrt(num_lipid[0] * BOXLX / BOXLY) + 0.5), ngridy = int(ngridx * BOXLY / BOXLX + 0.5);
  const double dx = BOXLX / ngridx, dy = BOXLY / ngridy;
  int o = 0;
  for (int k = 0; k < 2; ++k) {
    int ix = 0, iy = 0;
    double r0[3];
    r0[2] = z[k];
    for (int i = 0; i < num_lipid[k]; ++i) {
      if (ix == ngridx - 1) {
        ix = 0;
        ++iy;
      }
      r0[0] = (ix + 0.05 * (prng[0].gen_open0_open1() - 0.5)) * dx;
      r0[1] = (iy + 0.05 * (prng[0].gen_open0_open1() - 0.5)) * dy;
      // Three head beads
      double dr[] = { 0.0, 0.0, (k != 0) ? -L0_BOND_LIPID : L0_BOND_LIPID };
      for (int dim = 0; dim < 3; ++dim) 
        pos[o][dim] = r0[dim];
      pos[o+1][0] = pos[o][0] - 0.5 * L0_BOND_LIPID;
      pos[o+1][1] = pos[o][1] - HALF_SQRT3 * L0_BOND_LIPID;
      pos[o+1][2] = pos[o][2];
      pos[o+2][0] = pos[o+1][0] + 0.5 * L0_BOND_LIPID;
      pos[o+2][1] = pos[o+1][1] - HALF_SQRT3 * L0_BOND_LIPID;
      pos[o+2][2] = pos[o+1][2];
      // Two tails
      int oo = o + 3;
      for (int k = 0; k < 2; ++k) {
        for (int dim = 0; dim < 3; ++dim)
          pos[oo][dim] = pos[o+k+1][dim] + dr[dim];
        for (int i = 1; i < N_BEADS_PER_TAIL; ++i) {
          for (int dim = 0; dim < 3; ++dim)
            pos[oo+i][dim] = pos[oo+i-1][dim] + dr[dim];
        }
        oo += N_BEADS_PER_TAIL;
      }
      o += N_BEADS_PER_LIPID;
      ++ix;
    }
  }
  // Add water beads
  const int Nw = N_SOLVENT / 2;
  for (int i = N_LIPID_BEADS; i < N_BEADS; ++i) {
    pos[i][0] = BOXLX * prng[0].gen_open0_open1();
    pos[i][1] = BOXLY * prng[0].gen_open0_open1();
    if (i < Nw)
      pos[i][2] = z[0] * prng[0].gen_open0_open1();
    else
      pos[i][2] = z[1] + (BOXLZ - z[1]) * prng[0].gen_open0_open1();
    apply_pbc(pos[i]);
  }
  double Impulse[] = { 0.0, 0.0, 0.0 };
  for (int i = 0; i < N_BEADS; ++i)
    for (int dim = 0; dim < 3; ++dim) {
      vel[i][dim] = prng[0].gen_gauss() - 0.5;
      Impulse[dim] += vel[i][dim];
    } // assign initial velocity with Maxwell-Boltzmann distribution
  Ekin = 0.0;
  for (int i = 0; i < N_BEADS; ++i)
    for (int dim = 0; dim < 3; ++dim) {
      vel[i][dim] -= Impulse[dim] / (double)N_BEADS;
      Ekin += vel[i][dim] * vel[i][dim];
    }
  rescale_Tkin();
}

void read_init_conf() {
  std::ifstream in("dpd.cpt");
  for (int n = 0; n < N_LIPID_BEADS; ++n) // solute
    in >> pos[n][0] >> pos[n][1] >> pos[n][2]
       >> vel[n][0] >> vel[n][1] >> vel[n][2]
       >> box[n][0] >> box[n][1] >> box[n][2];
  for (int n = N_LIPID_BEADS; n < N_BEADS; ++n) // solvent
    in >> pos[n][0] >> pos[n][1] >> pos[n][2]
       >> vel[n][0] >> vel[n][1] >> vel[n][2];
  in.close(); 
}

void init_syst() {
  for (int t = 0; t < N_THREADS; ++t) {
    prng[t].init(time(0) + t * 2017);
    for (int i = 0; i < 1048576; ++i)
      prng[t].gen_int32();
  }
  /*
     W     H     T
  W  25.0  30.0  75.0
  H  30.0  30.0  35.0
  T  75.0  35.0  10.0
  */
  dpd_f_par[0][0][0] = 25.0; // Water-Water
  dpd_f_par[1][1][0] = 30.0; // Head-Head
  dpd_f_par[2][2][0] = 10.0; // Tail-Tail
  dpd_f_par[0][1][0] = dpd_f_par[1][0][0] = 30.0;
  dpd_f_par[0][2][0] = dpd_f_par[2][0][0] = 75.0;
  dpd_f_par[1][2][0] = dpd_f_par[2][1][0] = 35.0;

  dpd_f_par[0][0][1] = 4.5;
  dpd_f_par[1][1][1] = 4.5;
  dpd_f_par[2][2][1] = 4.5;
  dpd_f_par[0][1][1] = dpd_f_par[1][0][1] = 4.5;
  dpd_f_par[0][2][1] = dpd_f_par[2][0][1] = 20.0;
  dpd_f_par[1][2][1] = dpd_f_par[2][1][1] = 9.0;

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
       dpd_f_par[i][j][2] = sqrt(2.0 * dpd_f_par[i][j][1] * k_BT / DPD_DT);
  int count_pair = 0;
  for (int n = 0; n < N_LIPIDS; ++n) {
    int offset = n * N_BEADS_PER_LIPID;
    bonded_pair[count_pair  ][0] = offset + 0; bonded_pair[count_pair  ][1] = offset + 1;
    bonded_pair[count_pair+1][0] = offset + 1; bonded_pair[count_pair+1][1] = offset + 2;
    bonded_pair[count_pair+2][0] = offset + 1; bonded_pair[count_pair+2][1] = offset + 3;
    bonded_pair[count_pair+3][0] = offset + 3; bonded_pair[count_pair+3][1] = offset + 4;
    bonded_pair[count_pair+4][0] = offset + 4; bonded_pair[count_pair+4][1] = offset + 5;
    bonded_pair[count_pair+5][0] = offset + 5; bonded_pair[count_pair+5][1] = offset + 6;
    bonded_pair[count_pair+6][0] = offset + 2; bonded_pair[count_pair+6][1] = offset + 7;
    bonded_pair[count_pair+7][0] = offset + 7; bonded_pair[count_pair+7][1] = offset + 8;
    bonded_pair[count_pair+8][0] = offset + 8; bonded_pair[count_pair+8][1] = offset + 9;
    bonded_pair[count_pair+9][0] = offset + 9; bonded_pair[count_pair+9][1] = offset + 10;
    count_pair += 10;
  }
  int count_triple = 0;
  for (int n = 0; n < N_LIPIDS; ++n) {
    int offset = n * N_BEADS_PER_LIPID;
    bonded_triple[count_triple  ][0] = offset + 1; bonded_triple[count_triple  ][1] = offset + 3; bonded_triple[count_triple  ][2] = offset + 4; 
    bonded_triple[count_triple+1][0] = offset + 2; bonded_triple[count_triple+1][1] = offset + 7; bonded_triple[count_triple+1][2] = offset + 8;
    bonded_triple[count_triple+2][0] = offset + 3; bonded_triple[count_triple+2][1] = offset + 4; bonded_triple[count_triple+2][2] = offset + 5; 
    bonded_triple[count_triple+3][0] = offset + 4; bonded_triple[count_triple+3][1] = offset + 5; bonded_triple[count_triple+3][2] = offset + 6;
    bonded_triple[count_triple+4][0] = offset + 7; bonded_triple[count_triple+4][1] = offset + 8; bonded_triple[count_triple+4][2] = offset + 9;
    bonded_triple[count_triple+5][0] = offset + 8; bonded_triple[count_triple+5][1] = offset + 9; bonded_triple[count_triple+5][2] = offset + 10;
    count_triple += 6;
  }
  if (RESTART)
    read_init_conf();
  else
    gen_init_conf();
  gen_cell_list();
  write_log();
  write_traj(0);
  compute_forces();
}

#ifdef TENSION
void compute_bilayer_tension(int OPT) {
  static int    K = 0;
  static double st[N_SAMPLES][9] = { 0.0 };
  switch (OPT) {
    case 1 :
      if (K < N_SAMPLES) {
        for (int t = 0; t < N_THREADS; ++t)
          for (int i = 0; i < 9; ++i)
            st[K][i] += stress_tensor[t][i] / (BOXLX * BOXLY);
        ++K;
      }
      break;
    case 2 :
      if (K > 0) {
        std::ofstream out("stress_tensor.out");
        out.precision(8);
        for (int k = 0; k < K; ++k)
          out << st[k][0] << ' ' << st[k][1] << ' ' << st[k][2] << ' '
              << st[k][3] << ' ' << st[k][4] << ' ' << st[k][5] << ' '
              << st[k][6] << ' ' << st[k][7] << ' ' << st[k][8] << ' ' << (st[k][0] + st[k][4]) * 0.5 - st[k][8] << std::endl;
        out.close();
       }
      break;
  }
}
#endif

void check_conserv() {
  double p[3] = { 0.0, 0.0, 0.0 }, f[3] = { 0.0, 0.0, 0.0 };
  for (int n = 0; n < N_BEADS; ++n)
    for (int dim = 0; dim < 3; ++dim) {
      p[dim] += vel[n][dim];
      f[dim] += fff[0][n][dim];
    }
  std::cout << p[0] << ' ' << p[1] << ' ' << p[2] << ' ' << f[0] << ' ' << f[1] << ' ' << f[2] << std::endl;
}

int main() {
  init_syst();
  cTimer timer;
  if (RELAX) {
    for (unsigned int k = 0; k < INIT_CYCLE; ++k) {
      integrate();
      if (0 == k % FREQ_RESCALE) rescale_Tkin();
      if (0 == k % FREQ_SORT) sort_solvent();
      if (0 == k % FREQ_LOG) {
        std::ofstream out("sys.log", std::ios::app);
        out.precision(3);
        out.setf(std::ios::fixed, std::ios::floatfield);
        out << std::setw(10) << std::right << timer.elapsed() << " s ---> " 
            << std::setw(7) << std::right << k << '/' << INIT_CYCLE << " relax cycles\n";
        out.close();
      }
    }
  }
  for (unsigned int k = 1; k <= PROD_CYCLE; ++k) {
    integrate();
    //if (0 == k % 100) check_conserv();
    if (0 == k % FREQ_SORT) sort_solvent();
    if (0 == k % FREQ_SAMPLE) {
      write_traj();
      #ifdef TENSION
      compute_bilayer_tension(1);
      #endif
    }
    if (0 == k % FREQ_CPT) save_check_point();
    if (0 == k % FREQ_LOG) {
      std::ofstream out("sys.log", std::ios::app);
      out.precision(3);
      out.setf(std::ios::fixed, std::ios::floatfield);
      out << std::setw(10) << std::right << timer.elapsed() << " s ---> " 
          << std::setw(7) << std::right << k << '/' << PROD_CYCLE << " prod cycles, Epot/NkT =  "
          << std::setw(6) << std::left << Epot[0][0] / N_BEADS << ", Tkin = " 
          << std::setw(6) << std::left << Ekin / (3 * (N_BEADS - 1) * k_BT) << '\n';
      out.close();
    }
  }
  save_check_point();
  #ifdef TENSION
  compute_bilayer_tension(2);
  #endif
  double t_total = timer.elapsed();
  std::ofstream out("sys.log", std::ios::app);
  out << "\n---===CPU time===---\n";
  out << "Total CPU time: " << t_total << " seconds\n";
  out << " computing DPD force: " << t_force << " ( " << t_force / t_total * 100.0 << " %)\n";
  out << " generating linked list: " << t_list << " ( " << t_list / t_total * 100.0 << " %)\n";
  out << " sorting solvent particles: " << t_sort << " ( " << t_sort / t_total * 100.0 << " %)\n";
  out << " integrating motion equation: " << t_integrate << " ( " << t_integrate / t_total * 100.0 << " %)\n";
  out.close();
  return 0;
}
