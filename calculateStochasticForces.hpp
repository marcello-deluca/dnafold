#ifndef DEFINED_CALCULATE_STOCHASTIC_FORCES
#define DEFINED_CALCULATE_STOCHASTIC_FORCES
#include "headers.hpp"
void calculateStochasticForces(std::vector<std::vector<double> > & randomComponent, double k_B, double temp, double gamma_trans, size_t n_part, std::default_random_engine  & generator, std::normal_distribution<double> & distribution, double gamma_rot, size_t t, double r_var){
  for (size_t i = 0; i < n_part; ++i){
    for (size_t j = 0; j < 3; ++j){
      randomComponent[i][j] = sqrt (2 * k_B * temp / gamma_trans) * distribution (generator);
    }

    for (size_t j = 3; j < 4; ++j){
      randomComponent[i][j] = sqrt(r_var) * distribution(generator);
      randomComponent[i][j+1] = 0;
      randomComponent[i][j+2] = 0;
    }
    //if (t % 20000 == 0){
    // std::cout << "Random force components on " << i << " at time " << t << ": {" << randomComponent[i][0] << ", " << randomComponent[i][1] << ", " << randomComponent[i][2] << ", " << randomComponent[i][3] << ", " << randomComponent[i][4] << ", "<< randomComponent[i][5] << "}.\n";
    //}
  }
}

#endif
