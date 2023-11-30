#include <math.h>
#include <vector>

#ifndef DEFINED_PROJ
#define DEFINED_PROJ
//projects one 3d vector onto another one.
//pv = projecting vector, tv = target vector
std::vector<double> proj(std::vector<double> pv, std::vector<double> tv){
  double pdt = pv[0]*tv[0] + pv[1]*tv[1] + pv[2]*tv[2];
  double mag_t_sq = pow(tv[0],2) + pow(tv[1],2) + pow(tv[2],2);
  double pdt_div_tsq = pdt/mag_t_sq;
  std::vector<double> answer(3,0);
  answer[0] = pdt_div_tsq * tv[0];
  answer[1] = pdt_div_tsq * tv[1];
  answer[2] = pdt_div_tsq * tv[2];
  return answer;
}
#endif
