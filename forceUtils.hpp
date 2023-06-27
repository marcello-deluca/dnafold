//std=c++17
#ifndef FORCE_UTILS
#define FORCE_UTILS
#include "position3D.hpp"
#include "simBox.hpp"
#include <cmath>
#include <vector>
#include <iostream>
#include <tuple>
#include <mutex>

void applyPBC(bool pbc, double & RXIJ, double & RYIJ, double &RZIJ, simBox<double> & sb){
  if (pbc){
    RXIJ -= sb.dimensions.x * std::round(RXIJ / sb.dimensions.x);
    RYIJ -= sb.dimensions.y * std::round(RYIJ / sb.dimensions.y);
    RZIJ -= sb.dimensions.z * std::round(RZIJ / sb.dimensions.z);
  }
}

double dist (double dx, double dy, double dz){
  return sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
}

void addPairForces(double idx1, double idx2, std::vector<std::vector<double> > & forces, double FORCE, double RXIJ, double RYIJ, double RZIJ, double DIST, std::vector<std::mutex> & mtx){
  double FCUT = 3000;
  if (FORCE > FCUT){
    //std::cout << "Force bouncing off of cutoff\n";
    FORCE = FCUT;
  } else if (FORCE < -FCUT){\
    //std::cout << "Force bouncing off of cutoff\n";
    FORCE = -FCUT;
  }
  double FDIST = FORCE / DIST;
  double FX = FDIST * RXIJ;
  double FY = FDIST * RYIJ;
  double FZ = FDIST * RZIJ;
  mtx[idx1].lock();
  forces[idx1][0] += FX;
  forces[idx1][1] += FY;
  forces[idx1][2] += FZ;
  mtx[idx1].unlock();
  mtx[idx2].lock();
  forces[idx2][0] -= FX;
  forces[idx2][1] -= FY;
  forces[idx2][2] -= FZ;
  mtx[idx2].unlock();
  
}

void distances(bool pbc, simBox<double> & sb, double RXI, double RYI, double RZI, double & RXIJ, double & RYIJ, double & RZIJ, double & DIST, std::vector<position3D<double> > & r, size_t idx){
  RXIJ = r[idx].x - RXI;
  RYIJ = r[idx].y - RYI;
  RZIJ = r[idx].z - RZI;
  if (pbc){
    applyPBC(pbc,RXIJ,RYIJ,RZIJ,sb);
  }
  DIST = dist(RXIJ,RYIJ,RZIJ);
}
#endif
