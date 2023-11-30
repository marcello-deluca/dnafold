#ifndef DEFINED_BUILD_INTERACTION_LIST
#define DEFINED_BUILD_INTERACTION_LIST

#include "Interaction.hpp"
#include "position3D.hpp"
#include "forceUtils.hpp"

void BuildInteractionList(std::vector<Interaction> & ListOut, std::vector<position3D<double> > & r, simBox<double> & sb, size_t n_part, bool pbc, double CutoffRadius){
  double DIST, RXI, RYI, RZI, RXIJ, RYIJ, RZIJ;
  ListOut.clear();
  for (size_t p = 0; p < n_part; ++p){
    VerletList vl(p);
    for (size_t q = p; q < n_part; ++q){
      RXI = r[p].x;
      RYI = r[p].y;
      RZI = r[p].z;
      distances(pbc, sb, RXI, RYI, RZI,RXIJ, RYIJ, RZIJ, DIST, r, q);
      if (DIST < CutoffRadius && p!=q){
	Interaction IX = Interaction(p,q);
        ListOut.push_back(IX);
      }
    }
  }
}



#endif
