#ifndef recordFrame_H
#define recordFrame_H

#include <vector>
#include "position3D.hpp"
#include <iostream>
#include <ostream>
#include <fstream>
#include "simBox.hpp"
#include <functional>
#include <numeric>
#include <csignal>
void recordFrame(std::string & OUT_TYPE_, std::ofstream & traj, std::vector<position3D<double> > & R, size_t n_part, size_t n_scaf, std::vector<int> & stapleNumbers_, size_t t_, simBox<double> & simbox, std::vector<size_t> & isBound){
  signal(SIGINT, SIG_IGN); //disable ctrl+C kill so frame can be completely written.
  double NewCentroidX(0), NewCentroidY(0), NewCentroidZ(0);
  bool CenterScaffold = false;
  bool recordICS = false;
  if (OUT_TYPE_=="LAMMPS"){
    traj << "ITEM: TIMESTEP\n";
    traj << t_ << "\n";
    traj << "ITEM: NUMBER OF ATOMS\n";
    traj << n_part << "\n";
    traj << "ITEM: BOX BOUNDS pp pp pp\n";
    traj << -simbox.dimensions.x/2 << " " << simbox.dimensions.x/2<<"\n";
    traj << -simbox.dimensions.y/2 << " " << simbox.dimensions.y/2<<"\n";
    traj << -simbox.dimensions.z/2 << " " << simbox.dimensions.z/2<<"\n";
    traj << "ITEM: ATOMS id type xs ys zs\n";
    size_t type = 0;
    std::vector<double> xpos(n_part, 0);
    std::vector<double> ypos(n_part, 0);
    std::vector<double> zpos(n_part, 0);
    for (size_t i = 0; i < n_part; i++) {     
      xpos[i] = R[i].x/simbox.dimensions.x+0.5;
      ypos[i] = R[i].y/simbox.dimensions.y+0.5;
      zpos[i] = R[i].z/simbox.dimensions.z+0.5;
      if (CenterScaffold){      
	// CORRECTION FOR SCAFFOLD OVERFLOWING THROUGH PERIODIC BOUNDARY
	if (i>0 && i<n_scaf){
	  double ScaffoldDeltaX = xpos[i]-xpos[i-1];
	  if (ScaffoldDeltaX>0.25){
	    xpos[i]-=1;
	  } else if (ScaffoldDeltaX <-.25){
	    xpos[i]+=1;
	  }
	  double ScaffoldDeltaY = ypos[i]-ypos[i-1];
	  if (ScaffoldDeltaY>0.25){
	    ypos[i]-=1;
	  } else if (ScaffoldDeltaY<-.25){
	    ypos[i]+=1;
	  }

	  double ScaffoldDeltaZ = zpos[i]-zpos[i-1];
	  if (ScaffoldDeltaZ>0.25){
	    zpos[i]-=1;
	  } else if (ScaffoldDeltaZ<-0.25){
	    zpos[i]+=1;
	  }
	}
	double centroidX = 0;
	double centroidY = 0;
	double centroidZ = 0;
	for (size_t i = 0; i < n_scaf;i++){
	  centroidX += xpos[i];
	  centroidY += ypos[i];
	  centroidZ += zpos[i];
	}
	centroidX /= n_scaf;
	centroidY /= n_scaf;
	centroidZ /= n_scaf;
	// STAPLES IF THEY ARE NOT LOCATED CORRECTLY
	for (size_t i = n_scaf; i < n_part; i++){
	  if (i>=n_scaf && i < n_part){
	    if (isBound[i]){
	      double ScafStapDeltaX = xpos[i]-xpos[isBound[i]];
	      double ScafStapDeltaY = ypos[i]-ypos[isBound[i]];
	      double ScafStapDeltaZ = zpos[i]-zpos[isBound[i]];
	      if (ScafStapDeltaX>0.25){
		xpos[i]-=1;
	      } else if (ScafStapDeltaX<-.25){
		xpos[i]+=1;
	      }
	      if (ScafStapDeltaY>.25){
		ypos[i]-=1;
	      } else if (ScafStapDeltaY<-.25){
		ypos[i]+=1;
	      }
	      if (ScafStapDeltaZ>.25){
		zpos[i]-=1;
	      } else if (ScafStapDeltaZ<-.25){
		zpos[i]+=1;
	      }
	    }
	  }
	  if (xpos[i]>centroidX+0.5){
	    xpos[i]-=1;
	  }
	  else if (xpos[i]<centroidX-0.5){
	    xpos[i]+=1;
	  }
	  if (ypos[i]>centroidY+0.5){
	    ypos[i]-=1;
	  }
	  else if (ypos[i]<centroidY-0.5){
	    ypos[i]+=1;
	  }
	  if (zpos[i]>centroidZ+0.5){
	    zpos[i]-=1;
	  } else if (zpos[i]<centroidZ-0.5){
	    zpos[i]+=1;
	  }
	}
	NewCentroidX = 0;
	NewCentroidY = 0;
	NewCentroidZ = 0;
	for (size_t i = 0; i < n_part; i++){
	  NewCentroidX += xpos[i];
	  NewCentroidY += ypos[i];
	  NewCentroidZ += zpos[i];
	}
	NewCentroidX /= n_part;
	NewCentroidY /= n_part;
	NewCentroidZ /= n_part;
      }
    }
    for (size_t i = 0 ; i < n_part; i++){
      if (i < n_scaf){
	type = 1;
      } else {
	type = stapleNumbers_[i]+2;
      }
      traj << i+1 << " " << type << " " << xpos[i] << " " << ypos[i] << " " << zpos[i];
      if (!recordICS){
	traj << "\n";
      } else {
	traj << " " << R[i].a[0] << " " << R[i].a[1] << " " << R[i].a[2] << " " << R[i].b[0] << " " << R[i].b[1] << " " << R[i].b[2] << " " << R[i].c[0] << " " << R[i].c[1] << " " << R[i].c[2] << "\n";
      }
    }
  }
  signal(SIGINT, SIG_DFL); //reenable ctrl+C kill
}

#endif
