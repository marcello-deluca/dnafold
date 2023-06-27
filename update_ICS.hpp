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

void update_ICS(std::vector<position3D<double> > & r, simBox<double> & sb, bool pbc){
  double A_OLD0,A_OLD1,A_OLD2;
  double R_FORW0, R_FORW1, R_FORW2, DIST, MI;
  std::vector<double> da(3,0);
  double BDA;
  
  for (size_t i = 0, sz = r.size(); i < sz-1; i++){
    size_t IF1 = i+1;
    A_OLD0 = r[i].a[0];
    A_OLD1 = r[i].a[1];
    A_OLD2 = r[i].a[2];
    
    R_FORW0 = r[IF1].a[0]-r[i].a[0];
    R_FORW1 = r[IF1].a[1]-r[i].a[1];
    R_FORW2 = r[IF1].a[2]-r[i].a[2];

    if (pbc){
      R_FORW0 -= sb.dimensions.x*std::round(R_FORW0/sb.dimensions.x);
      R_FORW1 -= sb.dimensions.y*std::round(R_FORW1/sb.dimensions.y);
      R_FORW2 -= sb.dimensions.z*std::round(R_FORW2/sb.dimensions.z);
    }
    
    DIST = sqrt(pow(R_FORW0,2)+pow(R_FORW1,2)+pow(R_FORW2,2));

    MI = 1/DIST;
    r[i].a[0]=MI*R_FORW0;
    r[i].a[1]=MI*R_FORW1;
    r[i].a[2]=MI*R_FORW2;

    da[0]=r[i].a[0]-A_OLD0;
    da[1]=r[i].a[1]-A_OLD1;
    da[2]=r[i].a[2]-A_OLD2;

    //comparing to old portion
    BDA = dotprod(r[i].b, da);
    r[i].b[0] -= BDA * A_OLD0;
    r[i].b[1] -= BDA * A_OLD1;
    r[i].b[2] -= BDA * A_OLD2;

    //comparing to new section
    BDA = dotprod(r[i].b,r[i].a);
    r[i].b[0] -= BDA * r[i].a[0];
    r[i].b[1] -= BDA * r[i].a[1];
    r[i].b[2] -= BDA * r[i].a[2];

    MI = magnitude(r[i].b);
    r[i].b[0] *= MI;
    r[i].b[1] *= MI;
    r[i].b[2] *= MI;

    r[i].c[0] = r[i].a[1]*r[i].b[2] - r[i].a[2]*r[i].b[1];
    r[i].c[1] = r[i].a[2]*r[i].b[0] - r[i].a[0]*r[i].b[2];
    r[i].c[2] = r[i].a[0]*r[i].b[1] - r[i].a[1]*r[i].b[0];
    
  }
}
#endif
