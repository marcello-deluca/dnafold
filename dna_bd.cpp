#include "forceUtils.hpp"
#include "parameters.hpp"
#include "importer_function.hpp"
#include "simInit.hpp"
#include "recordFrame.hpp"
#include "resetForces.hpp"
#include "resetStochasticForces.hpp"
#include "calculateStochasticForces.hpp"
#include "calculateBendingForcesParallel.hpp"
#include "integrateMotionParallel.hpp"
#include "makeConnectivityMatrix.hpp"
#include "calculateICSTorsionAndBending.hpp"
#include "calculatePairForcesNeighbor.hpp"
#include <thread>
#include <mutex>
#include "MakeVerletList.hpp"
#include "VerletList.hpp"
int main (int argc, char ** argv){
  bool prependTrajWithLabels = false;
  bool debug_time = false;
  srand(seed);
  if (argc!=3){
    std::cerr << "usage: dnaBD designfile.json trajectory_out.dat\n";
  }
  ///////////////////////////////////////////////////////////////////
  // CREATE OUTPUT FILES ///////////////////////////////////////////
  /////////////////////////////////////////////////////////////////
  //std::ofstream trajectory;
  std::string trajectory_file_name;
  trajectory_file_name.append(argv[2]);
  std::ofstream trajectory(trajectory_file_name);
  //trajectory.open (trajectory_file_name);
  //trajectory << "TESTING TRAJECTORY OUTPUT\n";
  if (prependTrajWithLabels){
    trajectory << "n box x y z ax ay az bx by bz cx cy cz\n";
  }
  //std::ofstream bindTimeOutput;
  std::string bind_time_name = trajectory_file_name;
  bind_time_name.append("_BINDTIMES.dat");
  std::ofstream btout(bind_time_name);
  if (!btout){
    std::cout << "\n File \"" << bind_time_name << "\" did not open." << std::endl;
    return 1;
  }
  std::cout << "outputting bind times to: " << bind_time_name << ".\n";
  //btout << "Bind Times: \n" << "Time   number bound\n";


  std::string bind_time_name2 = trajectory_file_name;
  bind_time_name2.append("_BINDTIMES_AVERAGABLE.dat");
  std::ofstream btout2(bind_time_name2);
  if (!btout2){
    //std::cout << "\n File " << std::quoted(outFileName) << " did not open" << std::endl;
    std::cout << "\n File \"" << bind_time_name2 << "\" did not open." << std::endl;
    return 1;
  }
  std::cout << "outputting averageable bind times to: " << bind_time_name2 << ".\n";
  //btout << "Bind Times: \n" << "Time   number bound\n";


  std::string FirstBindTime = trajectory_file_name;
  FirstBindTime.append("_FIRSTBINDTIMES.dat");
  std::ofstream fbtOut(FirstBindTime);
  if (!fbtOut){
    //std::cout << "\n File " << std::quoted(outFileName) << " did not open" << std::endl;
    std::cout << "\n File \"" << FirstBindTime << "\" did not open." << std::endl;
    return 1;
  }
  std::cout << "outputting first bind times for each scaf to: " << FirstBindTime << ".\n";
  fbtOut << "particleNumber firstEncounterTime\n";
  
  std::ofstream forceOutput;
  std::string forceOutputName;
  forceOutputName.append(trajectory_file_name);
  forceOutputName.append("_FORCES.dat");
  forceOutput.open(forceOutputName);

  std::ofstream torqueOutput;
  std::string torqueOutputName;
  torqueOutputName.append(trajectory_file_name);
  torqueOutputName.append("_TORQUES.dat");
  torqueOutput.open(torqueOutputName);
  
  std::ofstream angleOutput;
  std::string angleOutputName;
  angleOutputName.append(trajectory_file_name);
  angleOutputName.append("_ANGLES.dat");
  angleOutput.open(angleOutputName);  
  angleOutput << "N B APG X Y\n";


  ////////////////////////////////////////////////////////////////////////
  //END OF OUTPUT FILE GENERATION ///////////////////////////////////////
  //////////////////////////////////////////////////////////////////////

  
  ////////////////////////////////////////////////////////////////////////
  // FILE LOADING AND CONNECTIVITY MATRICES ETC /////////////////////////
  //////////////////////////////////////////////////////////////////////
  std::cout << "LOADING FILE " << argv[1] << ".\n";
  load_file(argv[1], SM, staple_connections, stapleNumbers, connectivity, n_bonds, n_staple_seqs, isCrossover, StrandNumber, n_nucleotides_per_bead);
  //std::cout << "Connectivity of 111 and 112: " << connectivity[111][112] << " " << connectivity[112][111] << ".\n";
  std::cout << "LOADED " << argv[1] << ".\n";

  
  // for (size_t i = 0; i < n_scaf; i++){
  //   SM[i][i+128]=1;
  // }

  //for (size_t i = n_scaf; i<n_part-1;i++){
  //   staple_connections[i][i+1]=1;
  // }
  
  for (size_t i = 0; i < n_part;i++){
    for (size_t j = 0; j < n_part; j++){
      if (SM[i][j]){
	belongsTo[i]=j;
	belongsTo[j]=i;
      }
    }
  }
  //Sheet Structure
  isCrossover[42]=0;
  isCrossover[41]=0;
  isCrossover[belongsTo[42]]=0;
  isCrossover[belongsTo[41]]=0;

  
  //4HB
  //isCrossover[223]=1;
  //isCrossover[belongsTo[223]]=1;

  makeConnectivityMatrix(n_scaf, connectivity);
  for(size_t i = 0; i < n_part; i++){
    for (size_t j = 0; j < n_part; j++){
      if (SM[i][j]){
	belongsTo[i]=j;
      }
    }
  }

  
  std::ofstream PerfectComplements;
  std::string PerfectComplementsName;
  PerfectComplementsName.append("complements.txt");
  PerfectComplements.open(PerfectComplementsName);
  PerfectComplements << "ScaffoldParticle PairedStapleParticle StrandNumber\n";
  for (size_t i = 0; i < n_scaf; i++){
    PerfectComplements << i << " " << belongsTo[i] << "\n";
  }
  PerfectComplements.close();





  
  
  std::string OUT_TYPE = "LAMMPS";
  /////////////////////////////////////////////////////////////////////
  // END OF FILE LOADING /////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////
  // PRINT BASIC SIMULATION PARAMS //////////////////////////////////
  //////////////////////////////////////////////////////////////////
  std::cout << "Gamma_trans = " << gamma_trans <<".\n";
  std::cout << "Viscosity = " << visc << ".\n";
  std::cout << "Temp = " << temp << ".\n";
  std::cout << "Gamma_rot = " << gamma_rot << ".\n";
  std::cout << "Using seed = " << seed << ".\n";
  std::cout << "Initializing particles...\n";
  ///////////////////////////////////////////////////////////////////
  // END OF PRINT BASIC SIMULATION PARAMS //////////////////////////
  /////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////
  // INIT /////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  if (!RESTART){
    simInit(r, n_scaf, simbox, beadDiameter, n_part, beadAxialSeparation, beadRadialSeparation, staple_connections, SM);
  } else {
    load_frame(lastconfname,xboxsize,yboxsize,zboxsize,NATOM,STEPNUMBER,r);
  }

  if (RESTART){
    simbox.dimensions.x = xboxsize;
    simbox.dimensions.y = yboxsize;
    simbox.dimensions.z = zboxsize;
  }
  //////////////////////////////////////////////////////////////////
  // END OF INIT //////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  
  /////////////////////////////////////////////////////////////////
  // GENERATE TOPOLOGY FILE //////////////////////////////////////
  ///////////////////////////////////////////////////////////////
  std::string topology_filename = trajectory_file_name;
  topology_filename.append("_TOPOLOGY");
  std::ofstream topf;
  std::string LAMMPS_fname = topology_filename;
  LAMMPS_fname.append(".dat");
  topf.open(LAMMPS_fname);
  topf << "# Atom and bond counts\n";
  topf << n_part << " atoms\n";
  topf << "166 bonds\n";
  topf << "0 angles\n";
  topf << "0 dihedrals\n\n";
  topf << n_staple_seqs + 1 << " atom types\n";
  topf << "3 bond types\n";
  topf << "0 angle types\n";
  topf << "0 dihedral types\n\n";  
  topf << "# periodic box information\n";
  topf << -simbox.dimensions.x/2 << " " << simbox.dimensions.x/2 << " xlo xhi\n";
  topf << -simbox.dimensions.y/2 << " " << simbox.dimensions.y/2 << " ylo yhi\n";
  topf << -simbox.dimensions.z/2 << " " << simbox.dimensions.z/2 << " zlo zhi\n\n";
  topf << "Atoms\n\n";
  for (size_t i = 0; i < n_part; i++){
    topf << "   " << i << "   0   ";
    if (i < n_scaf){
      topf << 1;
    } else {
      topf << 2+stapleNumbers[i];
    }
    topf << "   " << r[i].x << "   " << r[i].y << "   " << r[i].z << "\n";
  }

  topf << std::endl;
  topf << "Bonds\n\n"; 
  int bond_idx = 0;
  for (size_t i = 0; i < n_part; i++){
    for (size_t j = 0; j < n_part; j++){
      if (staple_connections[i][j]) { // TYPE 1: STAPLE BACKBONE
        bond_idx++;
        topf << "   " << bond_idx << "   2   " << i << "   " << j << "\n";
      }
      if (connectivity[i][j]) { // TYPE 3: SCAFFOLD BACKBONE
        bond_idx++;
        topf << "   " << bond_idx << "   3   " << i << "   " << j << "\n";
      }
    }
  }
  topf.close();


  std::cout << "printed topology file\n";
  // return EXIT_SUCCESS;
  /////////////////////////////////////////////////////////////////
  // END OF GENERATE TOPOLOGY FILE ///////////////////////////////
  ///////////////////////////////////////////////////////////////
  
  /////////////////////////////////////////////////////////////////
  //// ACTUAL SIMULATION STARTS HERE /////////////////////////////
  ///////////////////////////////////////////////////////////////
  std::cout << "Starting Simulation...\n";
  // SET CLOCKS
  auto start = std::chrono::high_resolution_clock::now();
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds> (stop-start);
  auto start_frame = std::chrono::high_resolution_clock::now();
  auto stop_frame = std::chrono::high_resolution_clock::now();
  auto duration_frame = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
  size_t t = 0;
  if (!RESTART){
    t = 0;
  } else {
    t = STEPNUMBER;
  }
  
  ////////////////////////////////////////////////////////////////
  // SIMULATION LOOPS ///////////////////////////////////////////
  //////////////////////////////////////////////////////////////
  while(t < lsim) {
    start_frame = std::chrono::high_resolution_clock::now();
    // RECORD TO OUTPUT TRAJECTORY 
    if (t % stepsPerFrame == 0){
      start = std::chrono::high_resolution_clock::now();
      recordFrame(OUT_TYPE, trajectory, r, n_part, n_scaf, stapleNumbers, t, simbox, isBound);
      stop = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::microseconds> (stop-start);
      // RECORD FRAME TIME
      if (1){
	//std::cout << "Time to record frame = " << duration.count() << ".\n";
      }
    }
    // SHRINKING PERIODIC BOX
    if (simbox.dimensions.x > CubicBoxSize*simbox_final_size_ratio){
      simbox.dimensions.x -= shrink_rate;
      simbox.dimensions.y -= shrink_rate;
      simbox.dimensions.z -= shrink_rate;
    }
    // RESET FORCES
    resetForces(forces, torques,  n_part);
    resetStochasticForces(randomComponent, n_part);
    // OBTAIN RANDOM FORCES FROM IMPLICIT SOLVENT    
    calculateStochasticForces(randomComponent, k_B, temp, gamma_trans, n_part, generator, distribution, gamma_rot, t, r_var);
    if (t%1000==0){
      btout2 << t << " " << n_bound_staples << "\n";
    }
    if (t%stepsPerNeighborListRefresh==0){
      MakeVerletList(neighbors, r, simbox, n_part,pbc,12);
    }
    //PAIR POTENTIALS - BOTH BACKBONE AND NONBONDED
    start = std::chrono::high_resolution_clock::now();
    std::vector<std::thread> ForceLoopsThreads;
    ForceLoopsThreads.push_back(std::thread(calculatePairForces,std::ref(forces),std::ref(r),std::ref(simbox),std::ref(pbc),std::ref(isCrossover),std::ref(beadAxialSeparation),std::ref(beadRadialSeparation),std::ref(epsilon),std::ref(sigma),std::ref(n_part),std::ref(n_scaf),std::ref(n_stap),std::ref(belongsTo),std::ref(isBound),std::ref(staple_connections), std::ref(r_cut),std::ref(CIRCULAR_SCAFFOLD),std::ref(btout),std::ref(t),std::ref(n_bound_staples),true,std::ref(StrandNumber),std::ref(fbtOut),std::ref(prevBound),std::ref(binding_energy_kcal_mol),std::ref(binding_distance_cutoff),std::ref(FORCED_BINDING),std::ref(mutexes), std::ref(neighbors)));
    ForceLoopsThreads.push_back(std::thread(calculateBendingForces,n_scaf,std::ref(isBound),std::ref(r), std::ref(forces), dsdna_lp, k_B*temp, beadAxialSeparation, std::ref(stapleNumbers), std::ref(isCrossover), std::ref(simbox), pbc, std::ref(mutexes)));
    for (auto & thread:  ForceLoopsThreads){
      while(!thread.joinable()){
      }
      thread.join();
    }
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds> (stop-start);
    if (t%1000==0 && debug_time){
      std::cout << "debug: time for first set of force calculations in microseconds: " << duration.count() << ".\n";
    }
    // RUNGE KUTTA: PRETEND TO INTEGRATE AND STORE THE RESULT IN R_PROPOSED.
    start = std::chrono::high_resolution_clock::now();
    integrateMotion(n_part, forces, torques,  dt, r, r_proposed, randomComponent, verbose, pbc, simbox, n_scaf, gamma_trans, t, isCrossover, isBound);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds> (stop-start);
    if (t%1000==0 && debug_time){
      std::cout << "debug: duration for integration of particles in time for half of a timestep (pre runge kutta): " << duration.count() << ".\n";
    }
    // RECALCULATE PAIR FORCES, TORSION, AND BENDING FORCES BASED ON R_PROPOSED
    start = std::chrono::high_resolution_clock::now();
    resetForces(forces_forw, torques_forw, n_part);
    std::vector<std::thread> RungeKuttaThreads;
    RungeKuttaThreads.push_back(std::thread(calculatePairForces,std::ref(forces),std::ref(r_proposed),std::ref(simbox),std::ref(pbc),std::ref(isCrossover),std::ref(beadAxialSeparation),std::ref(beadRadialSeparation),std::ref(epsilon),std::ref(sigma),std::ref(n_part),std::ref(n_scaf),std::ref(n_stap),std::ref(belongsTo),std::ref(isBound),std::ref(staple_connections), std::ref(r_cut),std::ref(CIRCULAR_SCAFFOLD),std::ref(btout),std::ref(t),std::ref(n_bound_staples),true,std::ref(StrandNumber),std::ref(fbtOut),std::ref(prevBound),std::ref(binding_energy_kcal_mol),std::ref(binding_distance_cutoff),std::ref(FORCED_BINDING),std::ref(mutexes), std::ref(neighbors)));
    RungeKuttaThreads.push_back(std::thread(calculateBendingForces,n_scaf,std::ref(isBound),std::ref(r_proposed), std::ref(forces), dsdna_lp, k_B*temp, beadAxialSeparation, std::ref(stapleNumbers), std::ref(isCrossover), std::ref(simbox), pbc, std::ref(mutexes)));
    for (auto & thread : RungeKuttaThreads){
      while(!thread.joinable()){
      }
      thread.join();
    }
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds> (stop-start);
    if (t%1000 == 0 && debug_time){
      std::cout << "debug: duration for repeat of force calculations per runge kutta: " << duration.count() << ".\n";
    }
    //THESE ARE JUST FOR CHECKING THAT THE SUM OF FORCES AND TORQUES IS ZERO.
    std::vector<double> totalf(3,0);
    std::vector<double> totaltq(3,0);
    //Average forces
    for (size_t i = 0; i < n_part; i++){
      for (size_t j = 0; j < 3; j++){
      	forces[i][j]  = (forces[i][j]  +  forces_forw[i][j]) / 2;
	totalf[j] += forces[i][j];
      }
      torques[i][0] = (torques[i][0] + torques_forw[i][0])/2;
      totaltq[0] += torques[i][0];
    }
    if (t%1000==0){
      std::cout << "t = " << t*dt*1E-9*3314 << " / " << lsim*dt*1E-9*3314 << ".\n";
      std::cout << "force frame " << t/10000 << ": " << totalf[0] << ", " << totalf[1] << ", " << totalf[2] << "\n";  // torque: " << totaltq[0] << ".\n";
    }
    // Finally, integrate the motion for real based on the averaged forces, and store the result back into the regular r
    integrateMotion(n_part, forces, torques,  dt, r, r, randomComponent, verbose, pbc, simbox, n_scaf, gamma_trans, t, isCrossover, isBound);
    stop_frame = std::chrono::high_resolution_clock::now();
    duration_frame = std::chrono::duration_cast<std::chrono::microseconds>(stop_frame-start_frame);
    if (t%1001==0){
      std::cout << "Complete Timestep for step number" << t << " = " << duration_frame.count() << " microseconds.\n";
      std::cout << "Current temperature: " << temp << " Kelvin.\n";
    }
    // WHITELAM
    if (temp > final_temp){
      temp -= annealing_rate;
      tm273 = temp - 273.15;
      visc = (3.245157366681122E-11 * pow(tm273,4) - 9.061289916246450E-09 * pow(tm273,3) + 9.845093457119180E-07 * pow(tm273,2) - 5.521101038709857E-05 * tm273 + 1.778370453327189E-03) * 1E3;
      gamma_trans = 4 * _PI * visc * (beadAxialSeparation+beadRadialSeparation)/4;
      Er_linker = 4 * _PI * visc * pow(beadAxialSeparation/2,2) * beadAxialSeparation;
    }
    t++;
  } // NOW LOOP BACK TO THE NEXT TIMESTEP
  trajectory.close();
  torqueOutput.close();
  forceOutput.close();
  btout.close();
  return EXIT_SUCCESS;
}