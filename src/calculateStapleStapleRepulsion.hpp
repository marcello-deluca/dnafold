#ifndef STAPLE_STAPLE_REPULSION
#define STAPLE_STAPLE_REPULSION
#include "position3D.hpp"
#include "forceUtils.hpp"
#include "potentials.hpp"

void calculateStapleStapleRepulsion(std::vector<position3D<double> > & r, std::vector<std::vector<double> > & forces, const double beadRadialSeparation, const double beadAxialSeparation, const size_t n_scaf, const size_t n_stap, const size_t n_part, double epsilon, double sigma){

  double RXI, RYI, RZI;
  double RXIJ, RYIJ, RZIJ;
  double DIST = 0;
  double FORCE = 0;
  
  for (size_t i = n_scaf; i < n_part; i++){
    RXI = r[i].x;
    RYI = r[i].y;
    RZI = r[i].z;
    for (size_t j = i; j < n_part; j++){
      RXIJ = r[j].x - RXI;
      RYIJ = r[j].y - RYI;
      RZIJ = r[j].z - RZI;
      DIST = dist(RXIJ, RYIJ, RZIJ);
    }
    FORCE = wca(epsilon/5, sigma, DIST);
  }
}
#endif
