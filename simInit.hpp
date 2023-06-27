#ifndef DEFINED_SIMINIT
#define DEFINED_SIMINIT
#include "headers.hpp"
#include "noOverlapsAtInit.hpp"

void simInit(std::vector<position3D<double> > & r, const size_t & n_scaf, simBox<double> & simbox, const double beadDiameter, const size_t n_part, const double beadAxialSeparation,  const double beadRadialSeparation, std::vector<std::vector<int> > &staple_connections, std::vector<std::vector<int> > & SM){
  std::cout << "Initializing scaffold beads\n";

  bool initRandom = true;

  if (initRandom){
    r[0].x=0;
    r[0].y=0;
    r[0].z=0;

    for (size_t i = 1; i < n_scaf; i++){
      double dirx = (double) (rand() % 100) / 100 - .5; //random number between -.5 and .5
      double diry = (double) (rand() % 100) / 100 - .5; //random number between -.5 and .5
      double dirz = (double) (rand() % 100) / 100 - .5; //random number between -.5 and .5

      //std::cout << "Direction of travel: " << dirx << ", " << diry << ", " << dirz << ".\n";
    
      double mag = sqrt(pow(dirx,2) + pow(diry,2) + pow(dirz,2));
    
      // Make it so that the magnitude of init direction = 1.
      dirx /= mag;
      diry /= mag;
      dirz /= mag;


      //std::cout << "Direction of travel: " << dirx << ", " << diry << ", " << dirz << ".\n";

      r[i].x = r[i-1].x+ dirx * beadAxialSeparation;
      r[i].y = r[i-1].y+ diry * beadAxialSeparation;
      r[i].z = r[i-1].z+ dirz * beadAxialSeparation;
    }
  
  } else {

    for (size_t i = 0; i < n_scaf; i++){
      r[i].x = -simbox.dimensions.x / 4 + i * beadAxialSeparation + 6;
      r[i].y = 0;
      r[i].z = 0;
    }
    
  }
  
  for (size_t i = n_scaf; i < n_part; i++){
    do{
        r[i].x = -simbox.dimensions.x/2 + ((double)rand())/RAND_MAX * simbox.dimensions.x;
        r[i].y = -simbox.dimensions.y/2 + ((double)rand())/RAND_MAX * simbox.dimensions.y;
        r[i].z = -simbox.dimensions.z/2 + ((double)rand())/RAND_MAX * simbox.dimensions.z;
    }
    while (!noOverlapsAtInit(i, 2*beadRadialSeparation, r));
  }

  
  for (size_t i = n_scaf; i < n_part; ++i){
    for (size_t j = i; j < n_part; ++j){
      if (staple_connections[j][i]){
	//std::cout << "there is a staple connection at " << j << ", " << i << ".\n";
	r[j].x = r[i].x;
	r[j].y = r[i].y + beadAxialSeparation-0.1;
	r[j].z = r[i].z; 
      }
    }
  }


  for (size_t i = 0; i < n_scaf; i++){
    r[i].setA(1,0,0);
    r[i].setB(0,1,0);
    r[i].setC(0,0,1);
  }
  
  for (size_t i = 0; i < n_part;++i){
    for (size_t j = 0; j < n_part;++j){
      if (SM[i][j]){
        SM[j][i]=SM[i][j];
      }
      if (staple_connections[i][j]){
	staple_connections[j][i] = staple_connections[i][j];
      }
    }
  }
  std::cout << "Initialized particles\n";
}
#endif
