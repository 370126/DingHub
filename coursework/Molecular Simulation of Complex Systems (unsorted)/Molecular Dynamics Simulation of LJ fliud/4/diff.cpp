/**********************************************************
 NVE Molecular Dynamics simulations of Lennard-Jones fluid
 Calculation of diffusion coefficient
 written by HU, Jinglei 2018.04.09
**********************************************************/
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <cmath>

#define SQR(x) ((x) * (x))

#define n_samples 20000 
#define n_particles 500
#define delt_t 0.01     // time interval between two samples

double pos[n_samples][n_particles][3], vel[n_samples][n_particles][3];

int main() {
  int nn;
  std::string c;
  std::ifstream in("traj.xyz");
  for (int t = 0; t < n_samples; ++t) {
    in >> nn;
    in >> c >> c;
    for (int n = 0; n < n_particles; ++n)
      in >> c >> pos[t][n][0] >> pos[t][n][1] >> pos[t][n][2] >> vel[t][n][0] >> vel[t][n][1] >> vel[t][n][2];
  } 
  in.close();
  std::ofstream out("diff.dat");
  for (int dt = 0; dt < n_samples; ++dt) {
    int cnt = 0;
    double msd[3] = { 0.0, 0.0, 0.0 }, acf[3] = { 0.0, 0.0, 0.0 };
    for (int t0 = 0; t0 < n_samples - dt; ++t0) {
      int t1 = t0 + dt;
      for (int n = 0; n < n_particles; ++n) {
        for (int dim = 0; dim < 3; ++dim) {
          msd[dim] += SQR(pos[t1][n][dim] - pos[t0][n][dim]);
          acf[dim] += vel[t1][n][dim] * vel[t0][n][dim];
        }
        cnt++;
      }
    }
    out << dt * delt_t << "  " << msd[0] / cnt << "  " << msd[1] / cnt << "  " << msd[2] / cnt
                       << "  " << acf[0] / cnt << "  " << acf[1] / cnt << "  " << acf[2] / cnt << "  " << cnt << '\n';
  }
  out.close();
  return 0;
}
