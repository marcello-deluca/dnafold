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
void calculateBondedForces(std::vector<std::vector<double> > & forces, std::vector<position3D<double> > & r, simBox<double> & simbox, const size_t n_part, const size_t n_scaf, const double r_cut, bool pbc, std::vector<int> & StrandNumber, std::vector<bool> & isCrossover, std::vector<size_t> & isBound, double beadAxialSeparation, std::vector<std::vector<int> > & StrandNumber){
  double k_stretch = 14 * .0138 * 300; // (7kT/nm) 
  double RXI = 0;
  double RYI = 0;
  double RZI = 0;
  double RXIJ = 0;
  double RYIJ = 0;
  double RZIJ = 0;
  double DIST = 0;
  double FORCE = 0;
  size_t nscm1 = n_scaf-1;
  size_t nscm2 = n_scaf-2;
  
  for (size_t i = 1; i < nscm1; i++){
    size_t ip1 = i+1;
    //CURRENT POSITION
    RXI = r[i].x;
    RYI = r[i].y;
    RZI = r[i].z;
    //FORWARD DISTANCE:
    distances(pbc, RXI,RYI,RZI,RXIJ,RYIJ,RZIJ,DIST,r,ip1);
    if (StrandNumber[i] != StrandNumber[ip1] && isBound[i] && isBound[ip1]){
      FORCE = simpleHarmonic(beadRadialSeparation, DIST, k_stretch);
    } else { //NOT CROSSOVER
      FORCE = simpleHarmonic(beadAxialSeparation, DIST, k_stretch);
    }
    addPairForces(i,i+1,forces,FORCE,RXIJ,RYIJ,RZIJ,DIST,mtx);  
    FORCE = simpleHarmonic(beadAxialSeparation, DIST, k_stretch);
    addPairForces(i,ip1,forces,FORCE,RXIJ,RYIJ,RZIJ,DIST);

    //CROSSOVER STABILIZATION
    if (i < nScm2){
      if (isBound[i] && isBound[i+1] && isBound[i+2]){ //make sure all 3 species are bound
	if (StrandNumber[i]!=StrandNumber[i+1] && StrandNumber[i+1]!=StrandNumber[i+2] && StrandNumber[i] != StrandNumber[i+2]){ //2 consecutive crossovers
	  eqdist = sqrt(2*pow(beadRadialSeparation,2));
	  distances(pbc,sb,RXI,RYI,RZI,RXIJ,RYIJ,RZIJ,DIST,r,i+2);
	  FORCE = simpleHarmonic(eqdist, DIST, k_stretch*4);
	  addPairForces(i,i+2,forces,FORCE,RXIJ,RYIJ,RZIJ,DIST, mtx);
	}
	// BEFORE CROSSOVER AND 2 BEADS AFTER
	else if ( isCrossover[i] && isCrossover[i+1] && StrandNumber[i+1]==StrandNumber[i+2] && StrandNumber[i]!=StrandNumber[i+1]){
	  eqdist = sqrt(pow(beadAxialSeparation,2)+pow(beadRadialSeparation,2));
	  distances(pbc,sb,RXI,RYI,RZI,RXIJ,RYIJ,RZIJ,DIST,r,i+2);
	  FORCE = simpleHarmonic(eqdist, DIST,k_stretch*4);
	  addPairForces(i,i+2,forces,FORCE,RXIJ,RYIJ,RZIJ,DIST, mtx);
	  if (t%1000==0){
	    if (verbose){
	      std::cout << "pair " <<  i  << ", " << i+2  <<  " crossover = "  <<  FORCE << " from dist " << DIST << ".\n"; 
	    }
	  }
	}
	// 2 BEADS BEFORE CROSSOVER AND 1 BEAD AFTER
	else if (isCrossover[i+1] && isCrossover[i+2] && StrandNumber[i]==StrandNumber[i+1] && StrandNumber[i+1]!=StrandNumber[i+2]){
	  eqdist = sqrt(pow(beadAxialSeparation,2)+pow(beadRadialSeparation,2));
	  distances(pbc,sb,RXI,RYI,RZI,RXIJ,RYIJ,RZIJ,DIST,r,i+2);
	  FORCE = simpleHarmonic(eqdist, DIST,k_stretch*4);
	  addPairForces(i,i+2,forces,FORCE,RXIJ,RYIJ,RZIJ,DIST, mtx);
	  if (t%1000==0 && verbose){
	    std::cout << "pair " <<  i  << ", " << i+1  <<  " ssdna forces = "  <<  FORCE << " from dist " << DIST << ".\n"; 
	  }
	}
      }
    }
  }

  for (size_t i = n_scaf; i < n_part-1; i++){
    if (staple_connections[i][i+1]){
      distances(pbc,sb,RXI,RYI,RZI,RXIJ,RYIJ,RZIJ,DIST,r,i+1);
      if (StrandNumber[i] != StrandNumber[i+1] && isBound[i] && isBound[i+1]){
	FORCE = simpleHarmonic(beadRadialSeparation,DIST,k_stretch);
      } else {
	if (isBound[i] && isBound[i+1]){
	  FORCE = simpleHarmonic(beadAxialSeparation,DIST,k_stretch);
	} else {
	  FORCE = simpleHarmonic(ssdna_dist, DIST, ssdna_k);
	}
      }
      addPairForces(i,i+1,forces,FORCE,RXIJ,RYIJ,RZIJ,DIST,mtx);
    }
  }

  for (size_t i = 1; i < n_scaf; ++i){
    IDX = belongsTo[i];
    if (IDX != -1){
      distances(pbc,sb,RXI,RYI,RZI,RXIJ,RYIJ,RZIJ,DIST,r,IDX);
      if (DIST > BIND_CUTOFF){
	if (isBound[i] && PRE_RUNGE_KUTTA){
	  n_bound_staples--;
	  btout << t << " " << n_bound_staples << " " << i << "\n";
	  std::cout << t << " " << n_bound_staples << " " << i << "\n";
	  isBound[i]=0;
	  isBound[IDX]=0;
	}
	if (FORCED_BINDING){
	  FORCE = FORCED_BINDING_F;
	} else {
	  FORCE = 0;
	}
      } else {
	if (!isBound[i] && PRE_RUNGE_KUTTA){
	  // record first encounter time
	  if (!prevBound[i]){
	    prevBound[i]=1;
	    fbtOut << i << " " << t << std::endl; 
	  } else {
	    // do nothing
	  }
	  n_bound_staples++;
	  btout << t << " " << n_bound_staples << " " << i <<  "\n";
	  std::cout << t << " " << n_bound_staples << " " << i << "\n";
	  isBound[i] = IDX;
	  isBound[IDX] = i;
	}	  
	FORCE = hybridizationPotential(DIST,bind_force,0,bind_dist);
      }
      addPairForces(i,IDX,forces,FORCE,RXIJ,RYIJ,RZIJ,DIST, mtx);
    }
  }
}
#endif
