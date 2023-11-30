#ifndef DEFINED_CALCULATE_TORSION
#define DEFINED_CALCULATE_TORSION
#include "headers.hpp"
#include "potentials.hpp"
#include "proj.hpp"
#include "magnitude.hpp"
#include "difference.hpp"
#include "dotprod.hpp"
#include "forceUtils.hpp"
#include "updateICS.hpp"
#include "updateAlphaBetaGamma.hpp"



double computeBeta(position3D<double> &p1, position3D<double> &p2){
  double ada;
  ada = dotprod(p1.a,p2.a);
  if (ada >  1){ada =  1;}
  if (ada < -1){ada = -1;}
  return acos(ada);
}

double computeAlpha(position3D<double> &p1, position3D<double> &p2, double Beta){
  double sb, f1, Ac, f2;
  sb = sin(Beta);
  if (Beta >= 1E-10){
    f1 = dotprod(p2.a, p1.b)/sb;
  } else {
    f1 = dotprod(p2.b, p1.b);
  }
  if (f1 >  1){ f1 =  1;}
  if (f1 < -1){ f1 = -1;}
  Ac = acos(f1);
  f2 = dotprod(p2.a, p1.c);
  if (f2 >= 0){
    return Ac;
  } else {
    return -Ac;
  } 
}
double computeGamma(position3D<double> &p1, position3D<double> &p2, double Alpha){
  double f1, ada, apg, f2;
  ada = dotprod(p1.a,p2.a);
  f1 = (dotprod(p1.b, p2.b) + dotprod(p1.c, p2.c)) / (1.0 + ada);
  if (f1 > 1){ f1 = 1; }
  if (f1 < -1){ f1 = -1; }
  apg = acos(f1);
  f2 = (dotprod(p1.c, p2.b) - dotprod(p1.b, p2.c)) / (1.0 + ada);
  if (f2 >= 0){
    return  apg - Alpha;
  } else {
    return -apg - Alpha;
  }
}

std::vector<double> computeAlphaBetaGamma(position3D<double> &p1, position3D<double> &p2){
  double Beta = computeBeta(p1,p2);
  double Alpha = computeAlpha(p1,p2,Beta);
  double Gamma = computeGamma(p1,p2,Alpha);
  std::vector<double> ans;
  ans.push_back(Alpha);
  ans.push_back(Beta);
  ans.push_back(Gamma);
  return ans;
}



void calculateTorsion(std::vector<std::vector<double> > & torques, std::vector<std::vector<double> > & forces,  std::vector<position3D<double> > & r, simBox<double> &simbox, const size_t n_scaf, bool pbc, std::vector<size_t> & isBound, double kbt, const double beadAxialSeparation, size_t t, std::ofstream & angs, size_t stepsperframe, std::vector<bool> & isCrossover, std::vector<int> & stapleNumbers, std::vector<int> & belongsTo){

  // init
  bool bt_verbose = false;
  double twist_per_bead = _PI/4;//_PI/6; // 8 bases of twist, in radians
  bool stabilize_crossovers = true;
  bool ics_manage = false;
  std::vector<double> Alpha(n_scaf-2,0), Beta(n_scaf-2,0), Gamma(n_scaf-2,0);
  std::vector<double> a_old(3), r_forw(3), length(n_scaf-1, 0), da(3);
  //size_t im1, ip1;
  double lo, g, s;//, c1, c2, g1, g2, s1, s2;
  std::vector<double> a_m(3), df(3), Stri(3), Strim1(3), Ai(3), Aim1(3), Bi(3), Bim1(3), Chi(3), Chim1(3), Zhi(3), Zhim1(3), t_temp(3), z(3);
  
  // params
  double twist_cutoff = _PI;
  double twist_cutoff_nick = _PI;  

  lo = beadAxialSeparation;
  //h = 0;//10 * kbt / pow(beadAxialSeparation,2);
  g = 0;// 50 * kbt / beadAxialSeparation;
  s = 431.70 *2 * kbt / beadAxialSeparation;
  
  //double FMAX = 300;

  // THIS PORTION IS USED TO CALCULATE ALPHA, BETA, AND GAMMA
  size_t nm1 = n_scaf-1;
  
  updateICS(r,simbox,pbc, n_scaf);

  if (ics_manage && t>1000){
    for (size_t i = 2; i < n_scaf-2; i++){
      if (isCrossover[i] && isCrossover[i+1] && isBound[i] && isBound[i+1]){
        //r[i].a = r[i-1].a;
	// r[i].b = r[i-1].b;
        //r[i].c = r[i-1].c;
        for (size_t j = 0; j < 3; j++){
          r[i+1].a[j] = - r[i-1].a[j];
	  r[i+1].b[j] = r[i-1].b[j];
	  r[i+1].c[j] = - r[i-1].c[j];
        }
      }
    }
  }

  updateAlphaBetaGamma(Alpha,Beta,Gamma,r, n_scaf);

  if (t %20000==0){
    for (size_t i = 0; i < n_scaf-2; i++){
      //std::cout << "A[" << i << "] = " << Alpha[i] << ".\n";
      //std::cout << "G[" << i << "] = " << Gamma[i] << ".\n";
      std::cout << "APG[" << i << "] = " << Alpha[i] + Gamma[i] << ".\n";
      std::cout << "B[" << i << "] = " << Beta[i] << ".\n";
    }
  }
  
  //for (size_t i = 0; i < n_scaf-2; i++){
  //  if (isnan(Alpha[i]) || isnan(Gamma[i])){
  //      Alpha[i]=0;
  //      Gamma[i]=0;
  //  }
  // }
  
  if (t%stepsperframe==0){
    for (size_t i = 0; i < n_scaf-2;i++){
     angs << i << " " << Beta[i] << "  " << Alpha[i] + Gamma[i] <<" .\n";
    }
  }

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////// NOW FOR FORCE AND TORQUE CALCULATION  //////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
  // HAVE ALL BETA VALUES NOW, CALCULATE BENDING FORCE
  //nm1 = n_scaf-1;
  //Rcut2 = pow(r_cut,2);

  // // //for stability, assume g1 = 1/l[0] if Beta < 1E-5.
  // // if (Beta[0] > 1E-10){
  // //   g1 = Beta[0] / (sin(Beta[0]) * length[0]);
  // // } else {
  // //   g1 = 1.00/length[0]; //limit of x/sin(x) for numerical stability
  // // }
  // // c1 = cos(Beta[0]);


  //std::cout << "g1: " << g1 << ". c1: " << c1 << ".\n";
  
  // // // FIRST PARTICLE A, B CALCULATIONS

  // // length[0] = magnitude(r[0].a);
  // // Ai[0] = g1 * (r[1].a[0] - c1 * r[0].a[0]);
  // // Ai[1] = g1 * (r[1].a[1] - c1 * r[0].a[1]);
  // // Ai[2] = g1 * (r[1].a[2] - c1 * r[0].a[2]);
  // // Bi[0] = 0;
  // // Bi[1] = 0;
  // // Bi[2] = 0;

  // // g1 = (Alpha[0]+Gamma[0]+twist_per_bead) * tan(0.5*Beta[0]) / length[0];
  // // c1 = cos(Alpha[0]);
  // // s1 = sin(Alpha[0]);

  
  // // Chi[0] = g1 * (c1 * r[0].c[0] - s1 * r[0].b[0]);
  // // Chi[1] = g1 * (c1 * r[0].c[1] - s1 * r[0].b[1]);
  // // Chi[2] = g1 * (c1 * r[0].c[2] - s1 * r[0].b[2]);

  // // Zhi[0] = 0;
  // // Zhi[1] = 0;
  // // Zhi[2] = 0;


  // // //std::cout << "Chi 0: " << Chi[0] << ", " << Chi[1] << ", " << Chi[2] << ".\n";
  
  // // Stri[0] = (length[0] - lo) * r[0].a[0];
  // // Stri[1] = (length[0] - lo) * r[0].a[1];
  // // Stri[2] = (length[0] - lo) * r[0].a[2];

  // // //std::cout << "Calculated Chi, Zhi, Stri, no segfault yet\n";
  // // //if (isBound[0] && isBound[1]){
  
  // // if (isBound[0] && isBound[1]){
  // //   // First particle forward bending force
  // //   df[0] = -g * (Ai[0] + Bi[0]);
  // //   df[1] = -g * (Ai[1] + Bi[1]);
  // //   df[2] = -g * (Ai[2] + Bi[2]);

  // //   //forces[0][0] += df[0];
  // //   //forces[0][1] += df[1];
  // //   //forces[0][2] += df[2];
    
  // //   if (bt_verbose && t%10000==0){
  // //     std::cout << "force on 0 due to bending: " << df[0] << ", " << df[1] << ", " << df[2] << ".\n";
  // //   }
	
  // //   // First particle forward twisting 
  // //   df[0] = s * (Chi[0] + Zhi[0]);
  // //   df[1] = s * (Chi[1] + Zhi[1]);
  // //   df[2] = s * (Chi[2] + Zhi[2]);



    
  // //   if (isnan(df[0])||isnan(df[1])||isnan(df[2])){

  // //     std::cout << "nan triggered from Chi or Zhi calculation\n";
  // //     std::cout << "df: " << df[0] << ", " << df[1] << ", " << df[2] << ".\n";
  // //     std::cout << "Chi: " << Chi[0] << ", " << Chi[1] << ", " << Chi[2] << ".\n";
  // //     std::cout << "Zhi: " << Zhi[0] << ", " << Zhi[1] << ", " << Zhi[2] << ".\n";
  // //     std::cout << "Chi components: g1 = " << g1 << ", c1 = " << c1 << "s1 = " << s1 << ".\n";
  // //     std::cout << "g1 is problem child, components: Alpha[0]=" << Alpha[0] << ", Gamma[0]=" << Gamma[0] << ", Beta[0]=" << Beta[0] << ", length[0]=" << length[0] << ".\n"; 
  // //   }
    
  // //   //forces[0][0] += df[0];
  // //   //forces[0][1] += df[1];
  // //   //forces[0][2] += df[2];

  // //   if (bt_verbose && t%1000==0){
  // //     std::cout << "force on 0 due to twisting: " << df[0] << ", " << df[1] << ", " << df[2] << ".\n";
  // //   }
  // // }
  
  // // // OTHER PARTICLES: i TO n_scaf-2 (last particle follows)
  // // for (size_t i = 1; i < n_scaf-2; i++){

  // //   length[i] = magnitude(r[i].a);

  // //   im1 = i - 1;
  // //   ip1 = i + 1;    
  // //   Aim1[0] = Ai[0];
  // //   Aim1[1] = Ai[1];
  // //   Aim1[2] = Ai[2];
  // //   Bim1[0] = Bi[0];
  // //   Bim1[1] = Bi[1];
  // //   Bim1[2] = Bi[2];
    
  // //   if ( Beta[i] >= 1.0E-10 ) {
  // //     g1 = Beta[i] / (sin(Beta[i]) * length[i]);
  // //   } else { 
  // //     g1 = 1.0 / length[i];
  // //   }
    
  // //   c1 = cos(Beta[i]);
  // //   Ai[0] = g1 * (r[ip1].a[0] - c1 * r[i].a[0]);
  // //   Ai[1] = g1 * (r[ip1].a[1] - c1 * r[i].a[1]);
  // //   Ai[2] = g1 * (r[ip1].a[2] - c1 * r[i].a[2]);

  // //   if (Beta[im1] >= 1.0E-10 ){
  // //     g2 = Beta[im1] / (sin(Beta[im1]) * length[i]);
  // //   } else {
  // //     g2 = 1.0 / length[i];
  // //   }

  // //   c2 = cos(Beta[im1]);
  // //   Bi[0] = g2 * (r[im1].a[0] - c2 * r[i].a[0]);
  // //   Bi[1] = g2 * (r[im1].a[1] - c2 * r[i].a[1]);
  // //   Bi[2] = g2 * (r[im1].a[2] - c2 * r[i].a[2]);

  // //   //matches
  // //   Chim1[0] = Chi[0];
  // //   Chim1[1] = Chi[1];
  // //   Chim1[2] = Chi[2];
  // //   //matches
  // //   Zhim1[0] = Zhi[0];
  // //   Zhim1[1] = Zhi[1];
  // //   Zhim1[2] = Zhi[2];

    
  // //   //matches
  // //   g1 = (Alpha[i] + Gamma[i] + twist_per_bead) * tan(0.5 * Beta[i]) / length[i];
  // //   c1 = cos(Alpha[i]);
  // //   s1 = sin(Alpha[i]);

  // //   //matches
  // //   g2 = (Alpha[im1] + Gamma[im1] + twist_per_bead) * tan(0.5 * Beta[im1]) / length[i];
  // //   c2 = cos(Gamma[im1]);
  // //   s2 = sin(Gamma[im1]);

  // //   //matches
  // //   Chi[0] = g1 * (c1 * r[i].c[0] - s1 * r[i].b[0]);
  // //   Chi[1] = g1 * (c1 * r[i].c[1] - s1 * r[i].b[1]);
  // //   Chi[2] = g1 * (c1 * r[i].c[2] - s1 * r[i].b[2]);


  // //   //std::cout << "Chi " << i << ": " << Chi[0] << ", " << Chi[1] << ", " << Chi[2] << ".\n";
    
  // //   //matches
  // //   Zhi[0] = g2 * (c2 * r[i].c[0] + s2 * r[i].b[0]);
  // //   Zhi[1] = g2 * (c2 * r[i].c[1] + s2 * r[i].b[1]);
  // //   Zhi[2] = g2 * (c2 * r[i].c[2] + s2 * r[i].b[2]);

  // //   //std::cout << "Zhi " << i << ": " << Zhi[0] << ", " << Zhi[1] << ", " << Zhi[2] << ".\n";
    
  // //   //std::cout << "Chi = [" << Chi[0] << ", " << Chi[1] << ", " << Chi[2] << ".\n";
  // //   //std::cout << "Zhi = [" << Zhi[0] << ", " << Zhi[1] << ", " << Zhi[2] << ".\n";
    
  // //   //matches
  // //   Strim1[0] = Stri[0];
  // //   Strim1[1] = Stri[1];
  // //   Strim1[2] = Stri[2];

  // //   //matches
  // //   Stri[0] = (length[i] - lo) * r[i].a[0];
  // //   Stri[1] = (length[i] - lo) * r[i].a[1];
  // //   Stri[2] = (length[i] - lo) * r[i].a[2];

  // //   if (!isCrossover[i-1] && !isCrossover[i] && isBound[i-1] && isBound[i]){
  // //     df[0] = - g * (-Aim1[0] - Bim1[0]);
  // //     df[1] = - g * (-Aim1[1] - Bim1[1]);
  // //     df[2] = - g * (-Aim1[2] - Bim1[2]);
  // //     if (bt_verbose && t%1000==0){
  // // 	std::cout << "force on 0 due to bending: " << df[0] << ", " << df[1] << ", " << df[2] << ".\n";
  // //     }
	
  // //     if (magnitude(df) > FMAX){
  // // 	//std::cout << "force exceeding fmax detected at " << i << " from g-based restoring force.\n";
  // // 	df[0]=FMAX;
  // // 	df[1]=FMAX;
  // // 	df[2]=FMAX;
  // //     }
      
  // //     //forces[i][0] += df[0];
  // //     //forces[i][1] += df[1];
  // //     //forces[i][2] += df[2];

  // //     if (bt_verbose && t%20000==0){
  // // 	std::cout << "force on " << i << " due to bending (backwards): " << df[0] << ", " << df[1] << ", " << df[2] << ".\n";
  // //     }

  // //     df[0] = - s * (Chim1[0] + Zhim1[0]);
  // //     df[1] = - s * (Chim1[1] + Zhim1[1]);
  // //     df[2] = - s * (Chim1[2] + Zhim1[2]);

  // //     if (magnitude(df)>FMAX){
  // // 	//std::cout << "force exceeding fmax detected at " << i << " from s-based restoring force.\n";
  // // 	df[0]=FMAX;
  // // 	df[1]=FMAX;
  // // 	df[2]=FMAX;
  // //     }
      
  // //     //forces[i][0] += df[0];
  // //     //forces[i][1] += df[1];
  // //     //forces[i][2] += df[2];
	
  // //     if (bt_verbose && t%20000==0){
  // // 	std::cout << "force on " << i << " due to twisting (backwards): " << df[0] << ", " << df[1] << ", " << df[2] << ".\n";
  // //     }
  // //   }

  // //   if (i<n_scaf-1){
  // //     if (!isCrossover[i+1] && !isCrossover[i] && isBound[i] && isBound[i+1]){
  // // 	// if ((!isCrossover[i+1] || (isCrossover[i+1] && !isCrossover[i])) && isBound[i] && isBound[i+1]){

  // // 	if (t%20000 == 0){
  // // 	  //std::cout << "applying non-crossover forwards bending forces on " << i << ".\n";
  // // 	}
  // // 	df[0] =  s * (Chi[0] + Zhi[0]);
  // // 	df[1] =  s * (Chi[1] + Zhi[1]);
  // // 	df[2] =  s * (Chi[2] + Zhi[2]);
  // // 	//forces[i][0] += df[0];
  // // 	//forces[i][1] += df[1];
  // // 	//forces[i][2] += df[2];

  // // 	if (bt_verbose && t%20000==0){
  // // 	  std::cout << "force on " << i << " due to twisting (forwards): " << df[0] << ", " << df[1] << ", " << df[2] << ".\n";
  // // 	} 

  // // 	df[0] = -g * (Ai[0] + Bi[0]);
  // // 	df[1] = -g * (Ai[1] + Bi[1]);
  // // 	df[2] = -g * (Ai[2] + Bi[2]);
  // // 	//forces[i][0] += df[0];
  // // 	//forces[i][1] += df[1];
  // // 	//forces[i][2] += df[2];

  // // 	if (bt_verbose && t%20000==0){
  // // 	  std::cout << "force on " << i << " due to bending (forwards): " << df[0] << ", " << df[1] << ", " << df[2] << ".\n";
  // // 	} 
  // //     }
  // //   }
  // // } // 1 TO N_SCAF-1 ENDS RIGHT HERE
  
  // // // LAST PARTICLE

  // // // THESE Ai AND Bi ARE A(i-1) AND B(i-1) IN THE CONTEXT OF OUR ALGORITHM

  // // // // DON'T ENFORCE ANYTHING ON LAST PARTICLE -- TOO MUCH STUFF GOING ON WITH IT
  // // if (isBound[n_scaf-1] && isBound[n_scaf-2]){
  // //    df[0] = g * (Ai[0] + Bi[0]);
  // //    df[1] = g * (Ai[1] + Bi[1]);
  // //    df[2] = g * (Ai[2] + Bi[2]);
  // //    //forces[n_scaf-1][0] += df[0] - s * (Chi[0] + Zhi[0]);
  // //    //forces[n_scaf-1][1] += df[1] - s * (Chi[1] + Zhi[1]);
  // //    //forces[n_scaf-1][2] += df[2] - s * (Chi[2] + Zhi[2]);
  // //    if (bt_verbose && t%20000==0){
  // //      std::cout << "force on " << n_scaf-1 << " due to bending: " << df[0] << ", " << df[1] << ", " << df[2] << ".\n";
  // //      std::cout << "force on " << n_scaf-1 << " due to twisting: " << -s*(Chi[0]+Zhi[0]) << ", " <<  -s*(Chi[1]+Zhi[1]) << ", " <<  -s*(Chi[2]+Zhi[2]) << ".\n";
  // //    }
  // // }







  
  // THIS IS THE PART THAT WE ARE ACTUALLY USING RIGHT NOW






  /////////////////////////////////////////////
  //first particle forward-torque/////////////
  ///////////////////////////////////////////
  
  if (isBound[0] && isBound[1]){
    if (stapleNumbers[isBound[0]]==stapleNumbers[isBound[1]]){
      twist_cutoff = twist_cutoff_nick;
    } else {
      twist_cutoff = twist_cutoff_nick;
    }
    if ( (Alpha[0]+Gamma[0]) > twist_cutoff){
      torques[0][0] += s * (twist_cutoff-twist_per_bead);
      //torques[1][0] -= s * (twist_cutoff);
    } else if ((Alpha[0] + Gamma[0]) < -twist_cutoff){
      torques[0][0] -= s * (twist_cutoff-twist_per_bead);
      //torques[1][0] += s * twist_cutoff;
    } else {
      torques[0][0] += s * (Alpha[0] + Gamma[0] - twist_per_bead);
      //torques[1][0] -= s * (Alpha[0] + Gamma[0] + twist_per_bead);
    }
    if (t % 20000==0 && bt_verbose){
      std::cout << "torque on 0 = " << s*(Alpha[0]+Gamma[0] - twist_per_bead) << ".\n";
    }
  }
  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////
  ///////////////////////////////////////////////////// LOOKS OK

  ////////////////////////////////////////////////////////
  // FORWARD TORQUE, i->nm2 /////////////////////////////
  //////////////////////////////////////////////////////
  for (size_t i = 1; i < n_scaf-2; i++){
    //forward-torque
    if (!(isCrossover[i+1] && isCrossover[i]) && isBound[i] && isBound[i+1]){
      if (stapleNumbers[isBound[i]]==stapleNumbers[isBound[i+1]]){
	twist_cutoff = twist_cutoff_nick;
      } else {
	twist_cutoff = twist_cutoff_nick;
      }
      if ((Alpha[i]+Gamma[i]) > twist_cutoff){
	torques[i][0]   += s * (twist_cutoff);
	//torques[i+1][0] -= s * (twist_cutoff);
      } else if ((Alpha[i] + Gamma[i]) < -twist_cutoff){
	torques[i][0]   -= s * twist_cutoff;
	//torques[i+1][0] += s * twist_cutoff;
      } else {
	double t_add = s*(Alpha[i] + Gamma[i] - twist_per_bead);
        torques[i][0]   += t_add;
	torques[i+1][0] -= t_add;
      }
    }
    //backward-torque
    // if (isBound[i] && isBound[i-1]){//(!(isCrossover[i-1] && isCrossover[i]) && isBound[i-1] && isBound[i]){
    //   if (stapleNumbers[isBound[i-1]]==stapleNumbers[isBound[i]]){
    // 	twist_cutoff=_PI;
    //   } else {
    // 	twist_cutoff = twist_cutoff_nick;
    //   }
      
    //   if ( (Alpha[i-1]+Gamma[i-1]) > twist_cutoff ){
    //     torques[i][0]   -= s * twist_cutoff;
    // 	//torques[i+1][0] += s * twist_cutoff;
    //   } else if ( (Alpha[i-1] + Gamma[i-1]) < -twist_cutoff){
    //     torques[i][0] += s * twist_cutoff;
    // 	//torques[i+1][0]
    //   } else {
    // 	double t_add = s*(Alpha[i-1] + Gamma[i-1] + twist_per_bead);
    //     torques[i][0] -= t_add;
    //   }
    //   if (t % 20000==0 && bt_verbose){
    // 	std::cout << "torque on " << i << " = " << -s*(Alpha[i-1]+Gamma[i-1] + twist_per_bead) << ".\n";
    //   }
    // }
  }

  nm1 = n_scaf-1;
  size_t nm2 = n_scaf-2;

  /////////////////////////////////////////////
  //last particle backward-torque ////////////
  ///////////////////////////////////////////
  if (isBound[nm1] && isBound[nm2]){
    if (stapleNumbers[isBound[nm1]]==stapleNumbers[isBound[nm2]]){ //NO NICK
      twist_cutoff = _PI;
    } else { // THERE IS A NICK
      twist_cutoff = twist_cutoff_nick;
    }
    if ( (Alpha[nm2]+Gamma[nm2]) > twist_cutoff){
      torques[nm1][0] -= s * twist_cutoff;
      //torques[nm2][0] += s * twist_cutoff;
    } else if ((Alpha[nm2] + Gamma[nm2]) < -twist_cutoff){
      torques[nm1][0] -= s * -twist_cutoff;
      //torques[nm2][0] += s * -twist_cutoff;
    } else {
      torques[nm1][0] -= s * (Alpha[nm2] + Gamma[nm2] - twist_per_bead);
      //torques[nm2][0] += s * (Alpha[n_scaf-3] + Gamma[n_scaf-3] + twist_per_bead);
    }
  }

  if (stabilize_crossovers) {

    double crossover_angle = 0;
    double k_crossover = 1 * kbt;
    std::vector<double> rip3(3,0);
    std::vector<double> bi(3,0);
    std::vector<double> ang(n_scaf-1,0);
    std::vector<double> projection(3,0);
    std::vector<double> proj_onto_a_plane(3,0);
    std::vector<double> path(3,0);
    std::vector<double> fdir(3,0);
    for (size_t i = 0; i < n_scaf-1; i++) {
      for (size_t j = 0; j < n_scaf-1; j++) {

	//TORQUES TRANSMITTED ACROS STAPLE CROSSOVERS
	int iSt = belongsTo[i];
	int jSt = belongsTo[j];
	
	if ( i!=41 && j>i && isCrossover[iSt] && isCrossover[jSt] && isBound[i] && isBound[j] && stapleNumbers[iSt]==stapleNumbers[jSt]){
	  if (t%20000==0){
	    std::cout << i << " and " << j << " currently have ICS bound by staples.\n";
	  }
	  //COMPUTE EULER ANGLES AND APPLY TORQUES
	  //NEED FLIPPED AXES
	  position3D<double> fakeICS = r[j];
	  fakeICS.setA(-fakeICS.a[0], -fakeICS.a[1], -fakeICS.a[2]);
	  fakeICS.setB(-fakeICS.b[0], -fakeICS.b[1], -fakeICS.b[2]);
	  if (t%20000==0){
	    //std::cout << "Here is fakeICS: " << fakeICS.a[0] << ", " << fakeICS.a[1] << ", " << fakeICS.a[2] << "; " << fakeICS.b[0] << ", " << fakeICS.b[1] << ", " << fakeICS.b[2] << "; " << fakeICS.c[0] << ", " << fakeICS.c[1] << ", " << fakeICS.c[2] << ".\n";

	    //std::cout << "And here is real ICS: " << r[j].a[0] << ", " << r[j].a[1] << ", " << r[j].a[2] << "; " <<  r[j].b[0] << ", " << r[j].b[1] << ", " << r[j].b[2] << "; " <<  r[j].c[0] << ", " << r[j].c[1] << ", " << r[j].c[2] << ".\n"; 
	  }
	  std::vector<double> EulerAngles = computeAlphaBetaGamma(r[i],fakeICS);

	  if (t%20000==0){
	    //std::cout << i << " and " << j << " form a complete crossover, their ICSes are: " << r[i].a[0] << ", " << r[i].a[1] << ", " << r[i].a[2] << "; " <<  r[i].b[0] << ", " << r[i].b[1] << ", " << r[i].b[2] << "; " <<  r[i].c[0] << ", " << r[i].c[1] << ", " << r[i].c[2] << ";;;;;;;; " <<   r[j].a[0] << ", " << r[j].a[1] << ", " << r[j].a[2] << "; " <<  r[j].b[0] << ", " << r[j].b[1] << ", " << r[j].b[2] << "; " <<  r[j].c[0] << ", " << r[j].c[1] << ", " << r[j].c[2] << ".\n";

	    std::cout << "Alpha + Gamma = " << EulerAngles[0] + EulerAngles[2] << ".\n";
	    //std::cout << "Euler angles are " << EulerAngles[0] << ", " << EulerAngles[1] << ", " << EulerAngles[2] << ".\n";
	    
	    //std::cout << "Torque applied is: " << s * (EulerAngles[0] + EulerAngles[2]) << ".\n";
	  }
	  
	  double t_add  = s * 5 *(EulerAngles[0]+EulerAngles[2]-crossover_angle);
	  torques[i][0] += t_add;
	  torques[j][0] += t_add;

	  // COMPUTE PROJECTION OF i ONTO j's B VECTOR AND HARMONICALLY CONSTRAIN IN LINE
	  rip3[0] = r[j].x-r[i].x;
	  rip3[1] = r[j].y-r[i].y;
          rip3[2] = r[j].z-r[i].z;
	  if (pbc){
     	    rip3[0] -= simbox.dimensions.x * std::round(rip3[0]/simbox.dimensions.x);
     	    rip3[1] -= simbox.dimensions.y * std::round(rip3[1]/simbox.dimensions.y);
     	    rip3[2] -= simbox.dimensions.z * std::round(rip3[2]/simbox.dimensions.z);
     	  }
	  bi[0] = r[i].b[0];
     	  bi[1] = r[i].b[1];
     	  bi[2] = r[i].b[2];
    	  projection = proj(rip3,bi);
	  path[0] = projection[0] - rip3[0];
    	  path[1] = projection[1] - rip3[1];
    	  path[2] = projection[2] - rip3[2];

	  if (dotprod(rip3,bi)>0){
	    
	    forces[j][0] += path[0] * k_crossover;
    	    forces[j][1] += path[1] * k_crossover;
    	    forces[j][2] += path[2] * k_crossover;
	    forces[i][0] -= path[0] * k_crossover;
    	    forces[i][1] -= path[1] * k_crossover;
    	    forces[i][2] -= path[2] * k_crossover;
	  } else {
	    if (t%20000==0){
	       std::cout << "particle " << i << " and " << i+3 << " are misaligned.\n";
	    }
	  }
	}



	// APPLIED 0 SO THIS LOOP NEVER RUNS
	if (0&&isCrossover[i+1] && isCrossover[i+2] && j==i+3 && isBound[i] && isBound[j]){
	  rip3[0] = r[j].x-r[i].x;
	  rip3[1] = r[j].y-r[i].y;
          rip3[2] = r[j].z-r[i].z;
	  if (pbc){
     	    rip3[0] -= simbox.dimensions.x * std::round(rip3[0]/simbox.dimensions.x);
     	    rip3[1] -= simbox.dimensions.y * std::round(rip3[1]/simbox.dimensions.y);
     	    rip3[2] -= simbox.dimensions.z * std::round(rip3[2]/simbox.dimensions.z);
     	  }

	  //   double rip3_mag = magnitude(rip3);
	  bi[0] = r[i].b[0];
     	  bi[1] = r[i].b[1];
     	  bi[2] = r[i].b[2];

	  //   r_dot_b = rip3[0]*bi[0]+rip3[1]*bi[1]+rip3[2]*bi[2];

	  //   if (r_dot_b > 1){ r_dot_b = 1;}
	  //   if (r_dot_b < -1){ r_dot_b = -1;}

	  //   ang[i] = acos(r_dot_b);

	  //   if ( ang[i] >= 1.0E-10 ) {
	  //     g1 = ang[i] / (sin(ang[i]) * rip3_mag);
	  //   } else { 
	  //     g1 = 1.0 / rip3_mag;
	  //   }

	  //   c1 = cos(ang[i]);
	  //   Ai[0] = g1 * (r[ip1].a[0] - c1 * r[i].a[0]);
	  //   Ai[1] = g1 * (r[ip1].a[1] - c1 * r[i].a[1]);
	  //   Ai[2] = g1 * (r[ip1].a[2] - c1 * r[i].a[2]);
	  //   df[0] = -k_crossover * (Ai[0]);
          //   df[1] = -k_crossover * (Ai[0]);
          //   df[2] = -k_crossover * (Ai[2]);
          //   forces[i][0] += df[0];
          //   forces[i][1] += df[1];
          //   forces[i][2] += df[2];
	  //   forces[j][0] -= df[0];
	  //   forces[j][1] -= df[1];
	  //   forces[j][2] -= df[2];
	  
	  //    //torques[i][0] += k_crossover*ang[i];
	  // }
	  
	  // end of beta calculation	  

	  
    	  projection = proj(rip3,bi);
	  path[0] = projection[0] - rip3[0];
    	  path[1] = projection[1] - rip3[1];
    	  path[2] = projection[2] - rip3[2];

	  if (dotprod(rip3,bi)>0){
	    forces[j][0] += path[0] * k_crossover;
    	    forces[j][1] += path[1] * k_crossover;
    	    forces[j][2] += path[2] * k_crossover;
	    forces[i][0] -= path[0] * k_crossover;
    	    forces[i][1] -= path[1] * k_crossover;
    	    forces[i][2] -= path[2] * k_crossover;
	  } else {
	    if (t%20000==0){
	       std::cout << "particle " << i << " and " << i+3 << " are misaligned.\n";
	    }
	  }
	  
	
      
    
	  //APPLIED ZERO SO LOOP NEVER RUNS
	if (0 && isBound[i] && isBound[j] && isCrossover[i+1] && isCrossover[i+2] && j==i+3) { //torque on particle i to align its B vector with particle i+3
    	    rip3[0] = r[j].x-r[i].x;
    	    rip3[1] = r[j].y-r[i].y;
    	    rip3[2] = r[j].z-r[i].z;
	    if (pbc) { //Minimum Image Convention
    	      rip3[0] -= simbox.dimensions.x * std::round(rip3[0]/simbox.dimensions.x);
    	      rip3[1] -= simbox.dimensions.y * std::round(rip3[1]/simbox.dimensions.y);
    	      rip3[2] -= simbox.dimensions.z * std::round(rip3[2]/simbox.dimensions.z);
    	    }
    	    // Project R(i->i+3) onto a(i)
    	    projection = proj(rip3, r[i].a);
    	    std::vector<double> proj_r_a_plane = difference(rip3, projection);

    	    // theta_{a->b} = acos((a dot b)/(mag a * mag b))
    	    double dotprod_rip3_Bi =  dotprod(proj_r_a_plane, r[i].b);
    	    double m1 = magnitude(proj_r_a_plane);
    	    double m2 = magnitude(r[i].b);
    	    double theta_a_b = acos((dotprod_rip3_Bi)/(m1*m2));

    	    // This should be the first of 3 applications of torque.
    	    //torques[i][0] -= s * (theta_a_b); 

    	    if (t%20000==0){
    	      std::cout << "Here's what we got for crossover torque: " << -s*theta_a_b << ".\n";
    	    }  
    	  } //end of isBound[i] && isBound[j] && 0 loop
        } //end of for j loop

      // if (t%20000==0){
      // 	    std::cout << "Here is how far particle " << i << " is from b: " << path[0] << ", " << path[1] << ", " << path[2] << ".\n";
      // 	    std::cout << "This is rip3 dot bi for particle " << i << " : " << dotprod(rip3,bi) << ".\n";
      // 	    std::cout << "Here is how much force is acting on particle " << i+3 << ": " << path[0]*k_crossover<< ", " << path[1]*k_crossover << ", " << path[2] * k_crossover << ".\n";
      // 	    std::cout << "Here is the resulting motion from this alignment force: " << path[0] * k_crossover * .005 / (4*_PI*1.77837/10*3) << path[1] * k_crossover * .005 / (4*_PI*1.77837/10*3) << path[2] * k_crossover * .005 / (4*_PI*1.77837/10*3) << ".\n";
      // }

      
      } //end of j loop
    } //end of i loop
  } //end of if(stabilize crossovers) loop
} //END OF FUNCTION

#endif
