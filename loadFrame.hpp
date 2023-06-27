#ifndef LOAD_FRAME_H
#define LOAD_FRAME_H
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include "position3D.hpp"

void load_frame(std::string filename,double & xboxsize, double & yboxsize, double & zboxsize, int & NATOM, size_t & STEPNUMBER, std::vector<position3D<double> > & R){
  std::ifstream infile(filename);
  std::string BURN;
  std::string NAME;
  double boxxmin,boxxmax, boxymin, boxymax, boxzmin, boxzmax;
  infile >> BURN >> BURN;
  infile >> STEPNUMBER;
  infile >> BURN >> BURN >> BURN >> BURN; //ITEM: NUMBER OF ATOMS
  infile >> NATOM;
  infile >> BURN >> BURN >> BURN >> BURN >> BURN >> BURN; //ITEM: BOX BOUNDS pp pp pp
  infile >> boxxmin >> boxxmax;
  infile >> boxymin >> boxymax;
  infile >> boxzmin >> boxzmax;
  xboxsize = boxxmax - boxxmin;
  yboxsize = boxymax - boxymin;
  zboxsize = boxzmax - boxzmin;
  infile >> BURN >> BURN >> BURN >> BURN >> BURN >> BURN >> BURN; //ITEM: ATOMS id type xs ys zs
  for (int i = 0; i < NATOM; ++i){
    double XI;
    double YI;
    double ZI;
    infile >> BURN >> BURN >> XI >> YI >> ZI >> BURN >> BURN >> BURN >> BURN >> BURN >> BURN >> BURN >> BURN >> BURN; // ATNUM ATTYPE X Y Z AX AY AZ BX BY BZ CX CY CZ
    R[i].set((XI-0.5)*xboxsize, (YI-0.5) * yboxsize, (ZI-0.5) * zboxsize);
  }
  std::cout << " x y z bounds: " << boxxmin << ", " << boxxmax << "; " << boxymin << ", " << boxymax << "; " << boxzmin << ", " << boxzmax << ".\n";

  for (int i = 0; i < NATOM; ++i){
  std::vector<double> POS = R[i].getCoords();
  std::cout << "positions[" << i << "]: " << POS[0] << " " << POS[1] << " " << POS[2] << "\n";
    
  }
}
#endif
