#ifndef DEFINED_INTEGRATE_MOTION
#define DEFINED_INTEGRATE_MOTION

#include "headers.hpp"

void integrateMotion(const size_t n_part, std::vector<std::vector<double> > & forces, std::vector<std::vector<double> > & torques, const double dt, std::vector<position3D<double> > & r_source, std::vector<position3D<double> > & r_dest, std::vector<std::vector<double> > & randomComponent, bool verbose, bool pbc, simBox<double> & simbox, const size_t n_scaf, double gamma_trans, size_t t, std::vector<bool> & isCrossover, std::vector<size_t> & isBound){
  std::vector<std::vector<double> > d_theta(n_scaf, std::vector<double> (3, 0));
  double rad_linker = 1.33;
  //double eta = 1.78E-3;//1.778370453327189E-3; // viscosity of water
  // double Er_linker = 4*_PI*eta*pow(rad_linker,3)*2.0;
  // double stochastic_rotation = 0;
  // double deterministic_rotation = 0;
  // double MAX_STOCHASTIC_ROTATION = _PI;
  //std::vector<std::vector<double> > an (n_scaf, std::vector<double> (3,0));
  //std::vector<std::vector<double> > bn (n_scaf, std::vector<double> (3,0));
  //std::vector<std::vector<double> > cn (n_scaf, std::vector<double> (3,0));
  // double wa, wb, wc, wa2, wb2, wc2, z2, z, czt, szt, Omczt;
  // double g1, g2, g3;
  // for (size_t i = 0; i < n_scaf; i++){
  //   stochastic_rotation = randomComponent[i][3] * sqrt(dt);
  //   deterministic_rotation = dt * torques[i][0] / Er_linker;
  //   if (stochastic_rotation > MAX_STOCHASTIC_ROTATION){
  //     stochastic_rotation = MAX_STOCHASTIC_ROTATION;
  //   } else if (stochastic_rotation < -MAX_STOCHASTIC_ROTATION){
  //     stochastic_rotation = -MAX_STOCHASTIC_ROTATION;
  //   }
  //   d_theta[i][0] = deterministic_rotation + stochastic_rotation;
  //   d_theta[i][1] = 0;//randomComponent[i][4];
  //   d_theta[i][2] = 0;//randomComponent[i][5];
  // }
  // for (size_t i = 0; i < n_scaf; i++){
  //   for (size_t j = 0; j < 3; j++){
  //    an[i][j] = r_source[i].a[j];
  //    bn[i][j] = r_source[i].b[j];
  //    cn[i][j] = r_source[i].c[j];
  //  }
  // }
  // for (size_t i = 0; i < n_scaf; i++){
  //   wa  = d_theta[i][0];
  //   wb  = d_theta[i][1];
  //   wc  = d_theta[i][2];    
  //   wa2 = pow(wa,2);
  //   wb2 = pow(wb,2);
  //   wc2 = pow(wc,2);
  //   z2 = wa2 + wb2 + wc2;
  //   z = sqrt(z2);
  //   czt = (double)cos((double) (z * dt));
  //   szt = (double)sin((double) (z * dt));
  //   Omczt = 1-czt;
  //   if (z2 > 0){ // I think this is the only case you will realistically encounter.
  //     // ROTATE A
  //     g1 = ((wb2 + wc2) * czt + wa2) / z2;
  //     g2 = wa * wb * Omczt / z2 + wc * szt / z;
  //     g3 = wa * wc * Omczt / z2 - wb * szt / z;
  //     an[i][0] = g1 * r_source[i].a[0] + g2 * r_source[i].b[0] + g3 * r_source[i].c[0];
  //     an[i][1] = g1 * r_source[i].a[1] + g2 * r_source[i].b[1] + g3 * r_source[i].c[1];
  //     an[i][2] = g1 * r_source[i].a[2] + g2 * r_source[i].b[2] + g3 * r_source[i].c[2];
  //     // ROTATE B
  //     g1 = wa * wb * Omczt / z2 - wc * szt / z;
  //     g2 = ((wa2 + wc2) * czt + wb2) / z2;
  //     g3 = wb * wc * Omczt / z2 + wa * szt / z;
  //     bn[i][0] = g1 * r_source[i].a[0] + g2 * r_source[i].b[0] + g3 * r_source[i].c[0];
  //     bn[i][1] = g1 * r_source[i].a[1] + g2 * r_source[i].b[1] + g3 * r_source[i].c[1];
  //     bn[i][2] = g1 * r_source[i].a[2] + g2 * r_source[i].b[2] + g3 * r_source[i].c[2];
  //     // ROTATE C
  //     g1 = wa * wb * Omczt / z2 + wb * szt / z;
  //     g2 = wb * wc * Omczt / z2 - wa * szt / z;
  //     g3 = ((wa2 + wb2) * czt + wc2) / z2;
  //     cn[i][0] = g1 * r_source[i].a[0] + g2 * r_source[i].b[0] + g3 * r_source[i].c[0];
  //     cn[i][1] = g1 * r_source[i].a[1] + g2 * r_source[i].b[1] + g3 * r_source[i].c[1];
  //     cn[i][2] = g1 * r_source[i].a[2] + g2 * r_source[i].b[2] + g3 * r_source[i].c[2];
  //   } else { // If by some freak accident there is exactly zero change between steps...
  //     an[i][0] = r_source[i].a[0];
  //     an[i][1] = r_source[i].a[1];
  //     an[i][2] = r_source[i].a[2]; 

  //     bn[i][0] = r_source[i].b[0];
  //     bn[i][1] = r_source[i].b[1];
  //     bn[i][2] = r_source[i].b[2];

  //     cn[i][0] = r_source[i].c[0];
  //     cn[i][1] = r_source[i].c[1];
  //     cn[i][2] = r_source[i].c[2];
  //   }
  // }
  //for (size_t i = 0; i < n_scaf; i++){
  //   for (size_t j = 0; j < 3; j++){
  //     r_dest[i].a[j] = an[i][j];
  //     r_dest[i].b[j] = bn[i][j];
  //    r_dest[i].c[j] = cn[i][j];
  //  }
  // }
  for (size_t i = 0; i < n_part ; ++i){
    double deterministic_dx = forces[i][0] * dt / gamma_trans;
    double deterministic_dy = forces[i][1] * dt / gamma_trans;
    double deterministic_dz = forces[i][2] * dt / gamma_trans;
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
