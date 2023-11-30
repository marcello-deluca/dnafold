#include "headers.hpp"
#ifndef DEFINED_NOOVERLAPSATINIT
#define DEFINED_NOOVERLAPSATINIT
bool noOverlapsAtInit(size_t i, double min_dist, std::vector<position3D<double> > & r){
   for (size_t j = 0; j < i; j++){
    double dx = r[i].x-r[j].x;
    double dy = r[i].y-r[j].y;
    double dz = r[i].z-r[j].z;
    double dist = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));
    if (dist < min_dist){
      return false;
    }
  }
  return true;
}

#endif
