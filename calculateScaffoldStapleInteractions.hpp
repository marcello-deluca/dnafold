#ifndef DEFINED_CALCULATE_SCAFFOLD_STAPLE_INTERACTIONS
#define DEFINED_CALCULATE_SCAFFOLD_STAPLE_INTERACTIONS
#include "headers.hpp"
#include "potentials.hpp"
#include "forceUtils.hpp"

void calculateScaffoldStapleInteractions(std::vector<std::vector<double> > & forces, std::vector<position3D<double> > & r, simBox<double> & simbox, const size_t n_part, const size_t n_scaf, const double RCUT, std::vector<std::vector<int> > & SM, std::vector<size_t> & isBound, bool pbc, double k_stretch, double sigma, double epsilon){
  double BIND_CUTOFF = 0;
  double DIST = 0;
  double FORCE = 0;
  double RXI = 0;
  double RYI = 0;
  double RZI = 0;
  double RXIJ = 0;
  double RYIJ = 0;
  double RZIJ = 0;
  //Iterate over scaffold beads
  for (size_t i = 0; i < n_scaf; i++){
    RXI = r[i].x;
    RYI = r[i].y;
    RZI = r[i].z;
    //Iterate over staple beads
    for (size_t j = n_scaf; j < n_part; j++){
      RXIJ = r[j].x - RXI;
      RYIJ = r[j].y - RYI;
      RZIJ = r[j].z - RZI;
      applyPBC(pbc, RXIJ, RYIJ, RZIJ, simbox);
      DIST = dist(RXIJ, RYIJ, RZIJ);      
      //CASE 1: MATCHED STAPLE MATRIX, CAPABLE OF BINDING
      if (SM[i][j]){ // CASE 1: MATCHED STAPLE MATRIX, CAPABLE OF BINDING
	if (DIST > BIND_CUTOFF){
	  FORCE = 0;
	} else {
	  isBound[i] = j;
	  isBound[j] = i;
	  FORCE = simpleHarmonic(0, DIST, k_stretch);

	  addPairForces(i,j,forces,FORCE, RXIJ, RYIJ, RZIJ, DIST);
	}
      } else { // CASE 2: UNMATCHED, EXC VOL
	if (DIST > RCUT){
	  FORCE = 0;
	} else {
	  FORCE = wca(epsilon/5, sigma, DIST);
	  addPairForces(i,j,forces,FORCE,RXIJ,RYIJ,RZIJ,DIST);
	}
      }
    }
  }
}
#endif
