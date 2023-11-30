#ifndef DEFINED_RESET_FORCES
#define DEFINED_RESET_FORCES
#include "headers.hpp"
void resetForces(std::vector<std::vector<double> > & forces, std::vector<std::vector<double> > & torques, size_t n_part){
  for (size_t i = 0; i < n_part; ++i){
    forces[i][0] = 0; //x component
    forces[i][1] = 0; //y component
    forces[i][2] = 0; //z component
    torques[i][0] = 0;
    torques[i][1] = 0;
    torques[i][2] = 0;
  }
  
}
#endif
