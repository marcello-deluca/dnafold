#include <math.h>
#include <vector>
#ifndef DEFINED_MAGNITUDE
#define DEFINED_MAGNITUDE
double magnitude (std::vector<double> in1){
  return sqrt(pow(in1[0],2) + pow(in1[1],2) + pow(in1[2],2));
}
#endif
