
#ifndef DEFINED_CALCULATE_DETERMINISTIC_FORCES
#define DEFINED_CALCULATE_DETERMINISTIC_FORCES
#include "headers.hpp"
#include "potentials.hpp"

void calculateStretchingForces(std::vector<std::vector<double> > & forces, std::vector<position3D<double> > & r,simBox<double> & simbox, size_t n_part, double r_cut, bool pbc, std::vector<std::vector<size_t> > connectivity, double beadRadialSeparation, double crossover_stiffness, double beadAxialSeparation, std::vector<std::vector<size_t> > SM, bool delay_binding, std::vector<size_t> isBound, std::vector<std::vector<size_t> > staple_connections, size_t n_scaf, double l_k, double gamma_trans,double simbox_final_size_ratio){

  //std::cout << "starting force computation for step.\n";
  double RXI;
  double RYI;
  double RZI;
  double RXIJ;
  double RYIJ;
  double RZIJ;
  double force;
  std::vector<double> FX;
  std::vector<double> FY;
  std::vector<double> FZ;
  for (size_t i = 0, size = n_part; i<size; ++i){
    RXI = r[i].x;
    RYI = r[i].y;
    RZI = r[i].z;
    for (size_t j = 0, size = n_part; j<size; ++j){ //first, calculate distances between particles in x,y,z
      RXIJ = r[j].x - RXI;
      RYIJ = r[j].y - RYI;
      RZIJ = r[j].z - RZI;
      if (pbc){ //employ minimum image convention if pbc turned on.
	RXIJ -= simbox.dimensions.x * std::round(RXIJ/simbox.dimensions.x);
	RYIJ -= simbox.dimensions.y * std::round(RYIJ/simbox.dimensions.y);
	RZIJ -= simbox.dimensions.z * std::round(RZIJ/simbox.dimensions.z);
      }
      if (RXIJ > r_cut || RYIJ > r_cut || RZIJ > r_cut || i==j){
	forces[i][0] += 0;
	forces[i][1] += 0;
	forces[i][2] += 0;
      }
      else { // This is when we jump into all of the DNA-specific stuff.
	force = 0; //reset force first.
	double dr = sqrt(RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ);
	 

	if (dr > r_cut){ //total distance further than r_cut
	  force = 0;
	}

	//else if (((i%6==5 && j==i+1)&&i<36)||((j%6==5&&i==j+1)&&j<36)){
	//  force = simpleHarmonic(beadRadialSeparation, dr, crossover_stiffness);
	//}

	//else if ((i>56&&i%4==1&&j==i+1)||(j>56&&j%4==1&&i==j+1)){
	//  force = simpleHarmonic(beadRadialSeparation, dr, crossover_stiffness);
	//}

	else if (connectivity[i][j]){ //direct connection on scaffold
	  force = simpleHarmonic(beadAxialSeparation, dr, 1/l_k);
	}      

	else if (SM[i][j]){ //closer than r_cut and a designated staple-scaffold pair
	  if (delay_binding && simbox.dimensions.x > (simbox_final_size_ratio*n_part*beadAxialSeparation+2*r_cut+0.1)){
	    force = 0;
	  }
	  else {
	    isBound[i] = j;
	    isBound[j] = i;
	    if (dr > 4){ //Reduced cutoff for bound staple-scaffold beads
	      force = 0;
	    } else { //Harmonic potential applied according to eq dist from shared beads.
	      force = simpleHarmonic((8-SM[i][j])/8*beadAxialSeparation, dr, 1/l_k);
	    }
	  }
	}

	else if (staple_connections[i][j]){ //designated a connection between two sections of a staple
	  force = simpleHarmonic(beadAxialSeparation, dr, 1/l_k);
	  //force = modifiedFENE(beadDiameter, dr);
	}

	else if (i >=n_scaf && j >= n_scaf){
	  force = modifiedRepulsive(beadRadialSeparation, dr, 1/l_k);
	}
	else if (i >= n_scaf || j >= n_scaf){ //staple-staple soft repulsive interaction only
	  //CHANGED THIS TO ZERO TESTING ON AUG 6 2020!!!
	  force = 0;//modifiedRepulsivePotential(beadAxialSeparation, dr);
	  //if (dr > beadDiameter){
	  //  force = 0;
	  // }
	  //else {
	  //  force = modifiedHarmonic(beadDiameter,dr, 1/l_k);
	  // }
	}

	else { //some other soft repulsive interaction the same now but might change.
	  force = modifiedRepulsive(beadRadialSeparation, dr, 1/l_k);
	}

	forces[i][0] += (double) force * RXIJ / dr / gamma_trans;
	forces[i][1] += (double) force * RYIJ / dr / gamma_trans;
	forces[i][2] += (double) force * RZIJ / dr / gamma_trans;  
      }
    }
    
    //std::cout << "force on " << i << ": {" << forces(i,0)<< ", "<<forces(i,1)<<", "<<forces(i,2)<<"}.\n";
  }

}
#endif
