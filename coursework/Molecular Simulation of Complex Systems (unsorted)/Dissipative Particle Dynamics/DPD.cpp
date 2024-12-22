/****************************************************************************************
  Dissipative Particle Dynamics simulations of water
  Copyright (c) 2017.06 Jinglei Hu < hujinglei _at_ nju.edu.cn
  Macro usage:
    _OPENMP: openmp parallelization
****************************************************************************************/
/* head_of_celler files */
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
#define N_THREADS                 2
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
#define INIT_CYCLE                5000
#define PROD_CYCLE                20000
#define FREQ_RESCALE              50
#define FREQ_SAMPLE               20
#define N_CPT                     10
#define FREQ_CPT                  (PROD_CYCLE / N_CPT)
#define N_SAMPLES                 (PROD_CYCLE / FREQ_SAMPLE)
#define FREQ_LOG                  1000

/* Global variables */
cPRNG<MT19937, NORMAL::ZIGGURAT> prng[N_THREADS];
const double BOXL[3] = { BOXLX, BOXLY, BOXLZ }, HALF_BOXL[3] = { BOXLX * 0.5, BOXLY * 0.5, BOXLZ * 0.5 };
bool         at_boundary[N_CELLS][13]; // true if the neighboring cell is at boundary, or false otherwise
int          cell_list[N_CELLS][13], head_of_cell[N_CELLS], cell_idx[N_BEADS], linked_list[N_BEADS], idx_of_occur[N_BEADS];
double       pos[N_BEADS][3], vel[N_BEADS][3], fff[N_THREADS][N_BEADS+6][3] = { 0.0 }, _pos[N_BEADS][3], _vel[N_BEADS][3], _fff[N_BEADS][3];
double       t_force = 0.0, t_list = 0.0, t_sort = 0.0, t_integrate = 0.0, Ekin = 0.0, Epot[N_THREADS][18] = { 0.0 };
double       pressure_tensor[N_THREADS][9+16] = { 0.0 };

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
  for (int n = 0; n < N_BEADS; ++n)
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
  out << "Water\n";
  out << "Temperature: " << k_BT << '\n';
  out << "Box: " << BOXLX << " * " << BOXLY << " * " << BOXLZ << '\n';
  out << "RHO: " << RHO << '\n';
  out << "Total Number of Beads: " << N_BEADS << '\n';
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
  out << "Running...\n";
  out.close();
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
  #pragma omp parallel for schedule(static, N_BEADS / N_THREADS) num_threads(N_THREADS)
  #endif
  for (int n = 0; n < N_BEADS; ++n)
    cell_idx[n] = int(pos[n][0] * INV_CELL_LX) + (int(pos[n][1] * INV_CELL_LY) + int(pos[n][2] * INV_CELL_LZ) * N_CELLS_Y) * N_CELLS_X;
  for (int n = 0; n < N_BEADS; ++n) {
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
  #pragma omp parallel for schedule(static, N_BEADS / N_THREADS)  num_threads(N_THREADS)
  #endif
  for (int n = 0; n < counts; ++n) { // counts = N_BEADS
    const int k = idx_of_occur[n];
    for (int dim = 0; dim < 3; ++dim) {
      _pos[n][dim] = pos[k][dim];
      _vel[n][dim] = vel[k][dim];
      _fff[n][dim] = fff[0][k][dim];
    }
  }
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static, N_BEADS / N_THREADS)  num_threads(N_THREADS)
  #endif
  for (int n = 0; n < N_BEADS; ++n) {
    const int offset = 0;
    for (int dim = 0; dim < 3; ++dim) {
      pos[n][dim] = _pos[n][dim];
      vel[n][dim] = _vel[n][dim];
      fff[0][n][dim] = _fff[n][dim];
    }
  }
  t_sort += timer.elapsed();
} // to ensure solvent particles close in simulation box are stored close in memory

void calc_dpd_force(int i, int j, int TID = 0, bool MinImg = true) {
  double rij[] = { pos[i][0] - pos[j][0], pos[i][1] - pos[j][1], pos[i][2] - pos[j][2] };
  if (MinImg) apply_min_img(rij);
  const double rsq = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];
  if (rsq < DPD_R0_SQ) {
    const double r = sqrt(rsq), w = 1.0 / r - 1.0;
    const double a_ij = 25.0, g_ij = 4.5, s_ij = sqrt(2.0 * k_BT * g_ij / DPD_DT);
    Epot[TID][0] += 0.5 * a_ij * (1.0 + rsq - r - r);
    const double vij[] = { vel[i][0] - vel[j][0], vel[i][1] - vel[j][1], vel[i][2] - vel[j][2] }, \
                 rv = rij[0] * vij[0] + rij[1] * vij[1] + rij[2] * vij[2], \
                 coeff = (a_ij - g_ij * w * rv + s_ij * prng[TID].gen_gauss()) * w, \
                 ff[] = { rij[0] * coeff, rij[1] * coeff, rij[2] * coeff };
    for (int dim = 0; dim < 3; ++dim) {
      fff[TID][i][dim] += ff[dim];
      fff[TID][j][dim] -= ff[dim];
    }
    for (int dim1 = 0; dim1 < 3; ++dim1)
      for (int dim2 = 0; dim2 < 3; ++dim2)
        pressure_tensor[TID][3*dim1+dim2] += ff[dim1] * rij[dim2];
  }
}

void compute_forces() {
  gen_linked_list();
  cTimer timer;
  for (int t = 0; t < N_THREADS; ++t) {
    Epot[t][0] = 0.0;
    for (int i = 0; i < 9; ++i)
      pressure_tensor[t][i] = 0.0;
    for (int k = 0; k < N_BEADS; ++k)
      for (int dim = 0; dim < 3; ++dim)
        fff[t][k][dim] = 0.0;
  }

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
        int j = head_of_cell[cellj];
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
  for (int i = 0; i < N_BEADS; ++i) {
    pos[i][0] = BOXLX * prng[0].gen_open0_open1();
    pos[i][1] = BOXLY * prng[0].gen_open0_open1();
    pos[i][2] = BOXLZ * prng[0].gen_open0_open1();
    apply_pbc(pos[i]);
  }
  double Impulse[] = { 0.0, 0.0, 0.0 };
  for (int i = 0; i < N_BEADS; ++i)
    for (int dim = 0; dim < 3; ++dim) {
      vel[i][dim] = prng[0].gen_gauss();
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
  for (int n = 0; n < N_BEADS; ++n)
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
  if (RESTART)
    read_init_conf();
  else
    gen_init_conf();
  gen_cell_list();
  write_log();
  compute_forces();
}

void compute_pressure(int OPT) {
  static int    K = 0;
  static double pt[N_SAMPLES][9] = { 0.0 };
  switch (OPT) {
    case 1 :
      if (K < N_SAMPLES) {
        for (int t = 0; t < N_THREADS; ++t)
          for (int i = 0; i < 9; ++i)
            pt[K][i] += pressure_tensor[t][i] / VOL;
        for (int n = 0; n < N_BEADS; ++n)
          for (int dim1 = 0; dim1 < 3; ++dim1)
            for (int dim2 = 0; dim2 < 3; ++dim2)
              pt[K][dim1*3+dim2] += vel[n][dim1] * vel[n][dim2] / VOL;
        ++K;
      }
      break;
    case 2 :
      if (K > 0) {
        std::ofstream out("pressure_tensor.out");
        out.precision(8);
        for (int k = 0; k < K; ++k)
          out << pt[k][0] << ' ' << pt[k][1] << ' ' << pt[k][2] << ' '
              << pt[k][3] << ' ' << pt[k][4] << ' ' << pt[k][5] << ' '
              << pt[k][6] << ' ' << pt[k][7] << ' ' << pt[k][8] << std::endl;
        out.close();
       }
      break;
  }
}

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
    if (0 == k % FREQ_SAMPLE) compute_pressure(1);
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
  compute_pressure(2);
  double t_total = timer.elapsed();
  std::ofstream out("sys.log", std::ios::app);
  out.precision(3);
  out << "\n---===CPU time===---\n";
  out << "Total CPU time: " << t_total << " seconds\n";
  out << " computing DPD force: " << t_force << " ( " << t_force / t_total * 100.0 << " %)\n";
  out << " generating linked list: " << t_list << " ( " << t_list / t_total * 100.0 << " %)\n";
  out << " sorting solvent particles: " << t_sort << " ( " << t_sort / t_total * 100.0 << " %)\n";
  out << " integrating motion equation: " << t_integrate << " ( " << t_integrate / t_total * 100.0 << " %)\n";
  out.close();
  return 0;
}

