
#ifndef DEFINED_POTENTIALS
#define DEFINED_POTENTIALS
#include <cmath>
double hybridizationPotential(double dr, double force, double eq_dist, double sever_dist){
  if (dr > sever_dist){
    return 0;
  } else {
    return force;
  }
}

double simpleHarmonic(double eq_dist, double dist, double k_sp){
  return k_sp * (dist-eq_dist);
}

double modifiedRepulsive(double eq_dist, double dist, double k_sp){
  if (dist > eq_dist){
    return 0;
  } else {
    return k_sp * (dist-eq_dist);
  }
}
double wca(double sigma, double epsilon, double dist){
  if (dist > sigma * pow(2, 1/6)){
    return 0;
  } else {
    double fout = 4*epsilon*(-12 * pow(sigma,12) / pow(dist,13) + 6 * pow(sigma,6)/pow(dist,7));
    return fout;
  }
}

double excvolwca(double sigma, double kbt_k_ex, double dist){
  if (dist > sigma * pow(2,1/6)){
    return 0;
  } else {
    return 4 * kbt_k_ex * (-12 * pow(sigma,12) / pow(dist,13) + 6 * pow(sigma,6)/pow(dist,7)); 
  }
}
#endif
