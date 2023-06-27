#ifndef DEFINED_INTEGRATE_MOTION
#define DEFINED_INTEGRATE_MOTION

#include "headers.hpp"

void integrateMotionRange(size_t idx1, size_t idx2, const size_t n_part, std::vector<std::vector<double> > & forces, std::vector<std::vector<double> > & torques, const double dt, std::vector<position3D<double> > & r_source, std::vector<position3D<double> > & r_dest, std::vector<std::vector<double> > & randomComponent, bool verbose, bool pbc, simBox<double> & simbox, const size_t n_scaf, double gamma_trans, size_t t, std::vector<bool> & isCrossover, std::vector<size_t> & isBound){
  std::vector<std::vector<double> > d_theta(n_scaf, std::vector<double> (3, 0));
  for (size_t i = 0; i < n_part ; ++i){
    double deterministic_dx = forces[i][0] * dt / gamma_trans;
    double deterministic_dy = forces[i][1] * dt / gamma_trans;
    double deterministic_dz = forces[i][2] * dt / gamma_trans;
    // TO TEST SUPERCOILING: DO NOT ALLOW FIRST, SECOND, SECOND TO LAST,  OR LAST BEAD TO MOVE IN Y OR Z
    r_dest[i].x = (double) (r_source[i].x + deterministic_dx) + (double) randomComponent[i][0] * sqrt(dt);
    r_dest[i].y = (double) (r_source[i].y + deterministic_dy) + (double) randomComponent[i][1] * sqrt(dt);
    r_dest[i].z = (double) (r_source[i].z + deterministic_dz) + (double) randomComponent[i][2] * sqrt(dt);
    
  if (pbc){
      r_dest[i].x -= simbox.dimensions.x * std::round( r_dest[i].x / simbox.dimensions.x);
      r_dest[i].y -= simbox.dimensions.y * std::round( r_dest[i].y / simbox.dimensions.y);
      r_dest[i].z -= simbox.dimensions.z * std::round( r_dest[i].z / simbox.dimensions.z);      
    }
  }
}

void integrateMotion(const size_t n_part, std::vector<std::vector<double> > & forces, std::vector<std::vector<double> > & torques, const double dt, std::vector<position3D<double> > & r_source, std::vector<position3D<double> > & r_dest, std::vector<std::vector<double> > & randomComponent, bool verbose, bool pbc, simBox<double> & simbox, const size_t n_scaf, double gamma_trans, size_t t, std::vector<bool> & isCrossover, std::vector<size_t> & isBound){
  for (size_t i = 0; i < n_part ; ++i){
    double deterministic_dx = forces[i][0] * dt / gamma_trans;
    double deterministic_dy = forces[i][1] * dt / gamma_trans;
    double deterministic_dz = forces[i][2] * dt / gamma_trans;
    // TO TEST SUPERCOILING: DO NOT ALLOW FIRST, SECOND, SECOND TO LAST,  OR LAST BEAD TO MOVE IN Y OR Z
    r_dest[i].x = (double) (r_source[i].x + deterministic_dx) + (double) randomComponent[i][0] * sqrt(dt);
    r_dest[i].y = (double) (r_source[i].y + deterministic_dy) + (double) randomComponent[i][1] * sqrt(dt);
    r_dest[i].z = (double) (r_source[i].z + deterministic_dz) + (double) randomComponent[i][2] * sqrt(dt);
    
  if (pbc){
      r_dest[i].x -= simbox.dimensions.x * std::round( r_dest[i].x / simbox.dimensions.x);
      r_dest[i].y -= simbox.dimensions.y * std::round( r_dest[i].y / simbox.dimensions.y);
      r_dest[i].z -= simbox.dimensions.z * std::round( r_dest[i].z / simbox.dimensions.z);      
    }
  }
}
#endif
