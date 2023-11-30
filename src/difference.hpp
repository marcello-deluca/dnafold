#ifndef DEFINED_DIFFERENCE
#define DEFINED_DIFFERENCE

#include <vector>
#include <iostream>

std::vector<double> difference(std::vector<double> in, std::vector<double> to_subtract){
  std::vector<double> out;
  if (in.size() != to_subtract.size()){
    std::cerr << "attempted to take difference between two vectors of different sizes\n"; 
  } else {
    for (size_t i = 0; i < in.size();i++){
      out.push_back(in[i]-to_subtract[i]);
    }
  }
  return out;
}

#endif
