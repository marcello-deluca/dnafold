//-std=c++11
#ifndef DEFINED_CALCULATE_BENDING_FORCES
#define DEFINED_CALCULATE_BENDING_FORCES
#include "headers.hpp"
#include "simBox.hpp"
#include <iostream>
#include <thread>
#include <mutex>

void calculateIndividualBendingForces(size_t idx1, size_t idx2, const size_t n_scaf, std::vector<size_t> & isBound, std::vector<position3D<double> > & r, std::vector<std::vector<double> > & forces, double dsdna_lp, double kbt, const double beadAxialSeparation, std::vector<int> & stapleNumbers, std::vector<bool> & isCrossover, simBox<double> & sb, bool pbc, std::vector<std::mutex> & mtx){
  double THETA0 = 0; //default angle
  double KTH = dsdna_lp * kbt / beadAxialSeparation;
  for (size_t i = idx1; i<=idx2; i++){
    if ((isCrossover[i] && !isCrossover[i+1] && !isCrossover[i+2]) || (!isCrossover[i] && !isCrossover[i+1] && !isCrossover[i+2]) || (!isCrossover[i] && !isCrossover[i+1] && isCrossover[i+2]) || (isCrossover[i] && !isCrossover[i+1] && isCrossover[i+2]) || (isCrossover[i] && !isCrossover[i+1] && isCrossover[i+2]) ){
      KTH = dsdna_lp * kbt / beadAxialSeparation;
      THETA0 = 0;
    } else if ((isCrossover[i] && isCrossover[i+1]) || (isCrossover[i+1] && isCrossover[i+2]) || (isCrossover[i] && isCrossover[i+1] && isCrossover[i+2])){
      KTH = dsdna_lp * 25 * kbt/beadAxialSeparation;
      THETA0 = _PI/2;
    } else {
      std::cout << "Unexpected case encountered\n";
      std::cout << "case " << i << "\n";
      if (isCrossover[i]){
	std::cout << "i is crossover\n";
      }
      if (isCrossover[i+1]){
	std::cout << "i+1 is crossover\n";
      }
      if (isCrossover[i+2]){
	std::cout << "i+2 is crossover\n";
      }
    }
    if (isBound[i] && isBound[i+1] && isBound[i+2]){
      size_t s1 = stapleNumbers[isBound[i]];
      size_t s2 = stapleNumbers[isBound[i+1]];
      size_t s3 = stapleNumbers[isBound[i+2]];
      //if (s1==s2 && s2==s3){cutoffVal = 3.15;} else {cutoffVal = 1.5;}
      std::vector<double> ATOM1 = r[i].getCoords();
      std::vector<double> ATOM2 = r[i+1].getCoords();
      std::vector<double> ATOM3 = r[i+2].getCoords();
      double RX21 = ATOM2[0] - ATOM1[0];
      double RY21 = ATOM2[1] - ATOM1[1];
      double RZ21 = ATOM2[2] - ATOM1[2];
      double RX32 = ATOM3[0] - ATOM2[0];
      double RY32 = ATOM3[1] - ATOM2[1];
      double RZ32 = ATOM3[2] - ATOM2[2];
      if (pbc){
        RX21 -= sb.dimensions.x * std::round(RX21 / sb.dimensions.x);
        RY21 -= sb.dimensions.y * std::round(RY21 / sb.dimensions.y);
        RZ21 -= sb.dimensions.z * std::round(RZ21 / sb.dimensions.z);
	RX32 -= sb.dimensions.x * std::round(RX32 / sb.dimensions.x);
        RY32 -= sb.dimensions.y * std::round(RY32 / sb.dimensions.y);
	RZ32 -= sb.dimensions.z * std::round(RZ32 / sb.dimensions.z);
      }
      double C22 = RX21*RX21 + RY21*RY21 + RZ21*RZ21;
      double C33 = RX32*RX32 + RY32*RY32 + RZ32*RZ32;
      double C32 = RX21*RX32 + RY21*RY32 + RZ21*RZ32;
      // double C3322SQ = pow(C33*C22, 2);

      //OLD
      double C3322SQ = pow(C33*C22, 0.5);

      double COSTHETA = C32/C3322SQ;
      double THETA = acos(COSTHETA);
      double th_eff = THETA-THETA0; 
      if (abs(THETA) > _PI){th_eff = abs(THETA-2*_PI);}
      if (1){//!(abs(th_eff) > cutoffVal)){
	double GRADRX1 = (C32 * RX21 / C22 - RX32) / C3322SQ;
	double GRADRY1 = (C32 * RY21 / C22 - RY32) / C3322SQ;
	double GRADRZ1 = (C32 * RZ21 / C22 - RZ32) / C3322SQ;
	double GRADRX2 = (RX32 - RX21 + C32 * RX32 / C33 - C32 * RX21 /C22) / C3322SQ;
	double GRADRY2 = (RY32 - RY21 + C32 * RY32 / C33 - C32 * RY21 /C22) / C3322SQ;
	double GRADRZ2 = (RZ32 - RZ21 + C32 * RZ32 / C33 - C32 * RZ21 /C22) / C3322SQ;
	double GRADRX3 = (RX21 - C32 * RX32 / C33) / C3322SQ;
	double GRADRY3 = (RY21 - C32 * RY32 / C33) / C3322SQ;
	double GRADRZ3 = (RZ21 - C32 * RZ32 / C33) / C3322SQ;
	double FFI = KTH * th_eff / sin (th_eff);
	if (abs(THETA)<0.001){FFI = 0;}
	else {
	  mtx[i].lock();
	  forces[i][0] += FFI*GRADRX1;
          forces[i][1] += FFI*GRADRY1;
          forces[i][2] += FFI*GRADRZ1;
	  mtx[i].unlock();
	  mtx[i+1].lock();
	  forces[i+1][0] += FFI*GRADRX2;
          forces[i+1][1] += FFI*GRADRY2;
          forces[i+1][2] += FFI*GRADRZ2;
	  mtx[i+1].unlock();
	  mtx[i+2].lock();
          forces[i+2][0] += FFI*GRADRX3;
          forces[i+2][1] += FFI*GRADRY3;
	  forces[i+2][2] += FFI*GRADRZ3;
	  mtx[i+2].unlock();
	}
      }
    } // end if (isBound) statement
  }  
}


void calculateBendingForces(const size_t n_scaf, std::vector<size_t> & isBound, std::vector<position3D<double> > & r, std::vector<std::vector<double> > & forces, double dsdna_lp, double kbt, const double beadAxialSeparation, std::vector<int> & stapleNumbers, std::vector<bool> & isCrossover, simBox<double> & sb, bool pbc, std::vector<std::mutex> & mtx){
  double THETA0 = 0; //default angle
  double cutoffVal = 3.15; // radians
  double KTH = dsdna_lp * kbt / beadAxialSeparation;
  size_t NTHREADS = 1;
  std::vector<std::thread> threads;
  //Decomposed by particle index
  size_t n_part_per_thread = ceil(n_scaf-2 / NTHREADS);
  size_t idx1;
  size_t idx2;
  for (size_t i = 0; i < NTHREADS; i++){
    //figure out which particles this thread will compute forces for
    idx1 = i * n_part_per_thread;
    idx2 = (i+1) * n_part_per_thread-1;
    if (idx2>n_scaf-2){
      idx2 = n_scaf-2;
    }
    threads.push_back(std::thread(calculateIndividualBendingForces,idx1, idx2, n_scaf, std::ref(isBound), std::ref(r), std::ref(forces), dsdna_lp, kbt, beadAxialSeparation, std::ref(stapleNumbers), std::ref(isCrossover), std::ref(sb),std::ref(pbc), std::ref(mtx)));
  }
  for (auto & thread  : threads){
    while (!thread.joinable()){
      // do nothing
    }
    thread.join();
  }
}
#endif
