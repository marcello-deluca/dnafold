#ifndef DEFINED_RESET_STOCHASTIC_FORCES
#define DEFINED_RESET_STOCHASTIC_FORCES
#include "headers.hpp"
void resetStochasticForces(std::vector<std::vector<double> > & randomComponent, size_t n_part){
  for (size_t i = 0; i < n_part; i++){
    for (size_t j = 0; j < 3; j++){
      randomComponent[i][j]=0;
    }
  }
}



#endif
