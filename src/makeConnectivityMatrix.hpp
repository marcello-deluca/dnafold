#include "headers.hpp"
#ifndef DEFINED_MAKECONNECTIVITYMATRIX
#define DEFINED_MAKECONNECTIVITYMATRIX
void makeConnectivityMatrix(const size_t n_scaf, std::vector<std::vector<int> > & connectivity){
  for (size_t i = 0; i < n_scaf; i++){
    for (size_t j = 0; j < n_scaf; j++){
      if (j==i+1||j==i-1){
	connectivity[i][j]=1;
      }
      else {
	connectivity[i][j]=0;
      }
    }
  }
}
#endif
