#ifndef SCAFFOLD_FORCES_BONDED
#define SCAFFOLD_FORCES_BONDED
//#include "headers.hpp"
#include "simBox.hpp"
#include "position3D.hpp"
#include "potentials.hpp"
#include "forceUtils.hpp"
#include <vector>

// This function calculates forces between all bonded scaffold pairs.
// Written independently because it is O(N) instead of O(N^2).
void calculateBondedScaffoldForces(std::vector<std::vector<double> > & forces, std::vector<position3D<double> > & r, simBox<double> & simbox, const size_t n_part, const size_t n_scaf, const double r_cut, bool pbc, std::vector<bool> & isCrossover, double beadAxialSeparation){

  // PRE-INITIALIZE TO ZERO TO GUARANTEE NO UNDEFINED BEHAVIOR

  double k_stretch = 14 * .0138 * 300; // (7kT/nm)
  
  double RXI = 0;
  double RYI = 0;
  double RZI = 0;
  double RXIJ = 0;
  double RYIJ = 0;
  double RZIJ = 0;
  double DIST = 0;
  double FORCE = 0;
  
  for (size_t i = 1; i < n_scaf-1; i++){
    //CURRENT POSITION
    RXI = r[i].x;
    RYI = r[i].y;
    RZI = r[i].z;
    
    // BACKWARD DISTANCE:
    size_t im1 = i-1;
    RXIJ = r[im1].x - RXI;
    RYIJ = r[im1].y - RYI;
    RZIJ = r[im1].z - RZI;
    applyPBC(pbc, RXIJ, RYIJ, RZIJ, simbox);
    DIST = dist(RXIJ,RYIJ,RZIJ);
    FORCE = simpleHarmonic(beadAxialSeparation, DIST, k_stretch);
    addPairForces(i,i-1,forces,FORCE,RXIJ, RYIJ, RZIJ, DIST);
    
    //FORWARD DISTANCE:
    size_t ip1 = i+1;
    RXIJ = r[ip1].x - RXI;
    RYIJ = r[ip1].y - RYI;
    RZIJ = r[ip1].z - RZI;
    applyPBC(pbc, RXIJ, RYIJ, RZIJ, simbox);
    DIST = dist(RXIJ,RYIJ,RZIJ);
    FORCE = simpleHarmonic(beadAxialSeparation, DIST, k_stretch);
    addPairForces(i,i+1,forces,FORCE,RXIJ,RYIJ,RZIJ,DIST);
  }
}
#endif
