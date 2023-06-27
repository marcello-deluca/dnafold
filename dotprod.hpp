#ifndef DOTPROD_H
#define DOTPROD_H
#include <vector>
#include <iostream>

// assume vectors are of size 3
double dotprod(std::vector<double> & in1, std::vector<double> & in2){
  return in1[0]*in2[0]+in1[1]*in2[1]+in1[2]*in2[2];
}
#endif
