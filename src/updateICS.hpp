#ifndef UPDATE_ICS
#define UPDATE_ICS

#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include "position3D.hpp"
#include "dotprod.hpp"
#include "magnitude.hpp"
#include "simBox.hpp"

void updateICS(std::vector<position3D<double> > & r, simBox<double> & sb, bool pbc, size_t n_scaf){
  bool debug = false;
  double A_OLD0,A_OLD1,A_OLD2;
  double R_FORW0, R_FORW1, R_FORW2, DIST, MI;
  std::vector<double> da(3,0);
  double BDA;
  size_t IF1;

  if (debug){
    std::cout << "New forward distance calculation\n";
  }
  
  for (size_t i = 0, sz = n_scaf-1; i < sz; i++){

    IF1 = i+1;
    A_OLD0 = r[i].a[0];
    A_OLD1 = r[i].a[1];
    A_OLD2 = r[i].a[2];
    R_FORW0 = r[IF1].x - r[i].x;
    R_FORW1 = r[IF1].y - r[i].y;
    R_FORW2 = r[IF1].z - r[i].z;

    if (debug){
      std::cout << "forward dists before pbc correction for " << i << ": " << R_FORW0 << ", " << R_FORW1 << ", " << R_FORW2 << ".\n";
    }
    
    if (pbc){
      R_FORW0 -= sb.dimensions.x * std::round(R_FORW0 / sb.dimensions.x);
      R_FORW1 -= sb.dimensions.y * std::round(R_FORW1 / sb.dimensions.y);
      R_FORW2 -= sb.dimensions.z * std::round(R_FORW2 / sb.dimensions.z);
    }

    if (debug){
      std::cout << "Here are forward distances for " << i << ": " << R_FORW0 << ", " << R_FORW1 << ", " << R_FORW2 << ".\n";
    }

    DIST = sqrt(pow(R_FORW0,2)+pow(R_FORW1,2)+pow(R_FORW2,2));

    if (debug){
      std::cout << "DIST = " << DIST << ".\n";
    }

    MI = 1/DIST;

    r[i].a[0] = MI * R_FORW0;
    r[i].a[1] = MI * R_FORW1;
    r[i].a[2] = MI * R_FORW2;

    da[0] = r[i].a[0] - A_OLD0;
    da[1] = r[i].a[1] - A_OLD1;
    da[2] = r[i].a[2] - A_OLD2;

    //dotprod between b and da
    BDA = dotprod(r[i].b, da);
    r[i].b[0] -= BDA * A_OLD0;
    r[i].b[1] -= BDA * A_OLD1;
    r[i].b[2] -= BDA * A_OLD2;

    //dotprod between b and new a
    BDA = dotprod(r[i].b,r[i].a);
    r[i].b[0] -= BDA * r[i].a[0];
    r[i].b[1] -= BDA * r[i].a[1];
    r[i].b[2] -= BDA * r[i].a[2];

    MI = 1/ magnitude(r[i].b);
    r[i].b[0] *= MI;
    r[i].b[1] *= MI;
    r[i].b[2] *= MI;

    //Now have a and b, take crossprod to get c
    r[i].c[0] = r[i].a[1]*r[i].b[2] - r[i].a[2]*r[i].b[1];
    r[i].c[1] = r[i].a[2]*r[i].b[0] - r[i].a[0]*r[i].b[2];
    r[i].c[2] = r[i].a[0]*r[i].b[1] - r[i].a[1]*r[i].b[0];
  }
}

#endif
