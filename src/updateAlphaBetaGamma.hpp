#ifndef UPDATE_ALPHA_BETA_GAMMA
#define UPDATE_ALPHA_BETA_GAMMA
#include <vector>
#include <cmath>
#include "position3D.hpp"
#include "dotprod.hpp"

void updateAlphaBetaGamma(std::vector<double> & Alpha, std::vector<double> & Beta, std::vector<double> & Gamma, std::vector<position3D<double> > & r, size_t n_scaf){
  double sb;
  double f1;
  double f2;
  double ada;

  // Need i= 1 thru n_scaf-2 inclusive (fortran) or 0 thru  n_scaf-2 exclusive (c++)
  for (size_t i = 0, sz = n_scaf-2; i<sz; ++i){
    //take acos(dotprod(a(i), a(i+1))) to get Beta angle
    ada =  dotprod(r[i].a,r[i+1].a);
    if (ada >  1){ada =  1;}
    if (ada < -1){ada = -1;}
    Beta[i] = acos(ada);
    sb = sin(Beta[i]);

    // take acos(dotprod(b(i),a(i+1))/sin(Beta)) to get Alpha.
    if (Beta[i] >= 1E-10){
      f1 = dotprod(r[i+1].a, r[i].b)/sb;
    } else {
      f1 = dotprod(r[i+1].b, r[i].b);
    }
    if (f1 >  1){ f1 =  1;}
    if (f1 < -1){ f1 = -1;}
    double Ac = acos(f1);

    
    f2 = dotprod(r[i+1].a, r[i].c);
    if (f2 >= 0){
      Alpha[i] = Ac;
    } else {
      Alpha[i] = -Ac;
    }

    //Now use this to get Alpha+Gamma (total twist)
    f1 = (dotprod(r[i].b, r[i+1].b) + dotprod(r[i].c, r[i+1].c)) / (1.0 + ada);
    if (f1 > 1){ f1 = 1; }
    if (f1 < -1){ f1 = -1; }
    double apg = acos(f1);
    f2 = (dotprod(r[i].c, r[i+1].b) - dotprod(r[i].b, r[i+1].c)) / (1.0 + ada);
    if (f2 >= 0){
      Gamma[i] =  apg - Alpha[i];
    } else {
      Gamma[i] = -apg - Alpha[i];
    }
  }
}
#endif
