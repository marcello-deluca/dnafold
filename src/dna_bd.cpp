#include "forceUtils.hpp"
#include "parameters.hpp"
#include "importer_function.hpp"
#include "simInit.hpp"
#include "recordFrame.hpp"
#include "setParameters.hpp"
#include "readInputFile.hpp"
#include "resetForces.hpp"
#include "resetStochasticForces.hpp"
#include "calculateStochasticForces.hpp"
#include "calculateBendingForces.hpp"
#include "integrateMotionParallel.hpp"
#include "makeConnectivityMatrix.hpp"
#include "calculateICSTorsionAndBending.hpp"
#include "calculatePairForcesNeighbor.hpp"
#include <thread>
#include <mutex>
#include "MakeVerletList.hpp"
#include "VerletList.hpp"
#include "parse_json.hpp"
int main (int argc, char ** argv){
  bool prependTrajWithLabels = false;
  bool debug_time = false;
  if (argc!=2){
    std::cerr << "usage: dnaBD inputFile\n";
  }
  ///////////////////////////////////////////////////////////////////
  // CREATE OUTPUT FILES ///////////////////////////////////////////
  /////////////////////////////////////////////////////////////////
  std::vector<std::pair<std::string, std::string> > inputs;
  readInputFile((std::string)argv[1], inputs);
  std::string inFile;
  std::string outFile;
  size_t NTHREADS = 1;
  int reportTimings = 0;
  // RESTART SPECS AND VARIABLES
  bool RESTART = false;
  std::string lastconfname = "last_conf.dat";
  double xboxsize, yboxsize, zboxsize;
  int NATOM;
  size_t STEPNUMBER = 0;
  setParameters(stepsPerFrame,lsim,dt,print_to_stdout_every,stepsPerNeighborListRefresh,temp,final_temp,annealing_rate,lastconfname,CIRCULAR_SCAFFOLD,FORCED_BINDING,NTHREADS,reportTimings,pbc, inFile,outFile, binding_energy_kcal_mol, CubicBoxSize, shrink_rate, simbox_final_size_ratio,inputs);
  simBox<double> simbox(CubicBoxSize, CubicBoxSize, CubicBoxSize);
  std::string trajectory_file_name;
  std::cout << "Outputting to files starting with " << outFile << "\n";
  trajectory_file_name.append(outFile);
  std::ofstream trajectory(trajectory_file_name);
  if (prependTrajWithLabels){
    trajectory << "n box x y z ax ay az bx by bz cx cy cz\n";
  }
  std::string bind_time_name = trajectory_file_name;
  bind_time_name.append("_BINDTIMES.dat");
  std::ofstream btout(bind_time_name);
  if (!btout){
    std::cout << "\n File \"" << bind_time_name << "\" did not open." << std::endl;
    return 1;
  }
  std::cout << "outputting bind times to: " << bind_time_name << ".\n";
  std::string bind_time_name2 = trajectory_file_name;
  bind_time_name2.append("_BINDTIMES_AVERAGABLE.dat");
  std::ofstream btout2(bind_time_name2);
  if (!btout2){
    //std::cout << "\n File " << std::quoted(outFileName) << " did not open" << std::endl;
    std::cout << "\n File \"" << bind_time_name2 << "\" did not open." << std::endl;
    return 1;
  }
  std::cout << "outputting averageable bind times to: " << bind_time_name2 << ".\n";
  std::string FirstBindTime = trajectory_file_name;
  FirstBindTime.append("_FIRSTBINDTIMES.dat");
  std::ofstream fbtOut(FirstBindTime);
  if (!fbtOut){
    std::cout << "\n File \"" << FirstBindTime << "\" did not open." << std::endl;
    return 1;
  }
  std::cout << "outputting first bind times for each scaf to: " << FirstBindTime << ".\n";
  fbtOut << "particleNumber firstEncounterTime\n";
  
  std::vector<std::vector<int> > scaffoldNucleotides;
  std::vector<std::vector<int> > stapleNucleotides;
  
  std::cout << "LOADING FILE " << inFile.data() << "\n";
  parse_json(&inFile[0], n_scaf, n_stap, scaffoldNucleotides, stapleNucleotides);

  size_t n_part = n_scaf + n_stap;
  std::cout << "Detected " << n_scaf << " scaffold beads and " << n_stap << " staple beads for a total of " << n_part << " beads in this design.\n";

  // INITIALIZE VECTORS
  std::vector<std::mutex> mutexes(n_part);
  std::vector<size_t> isBound(n_part, default_size_t);
  std::vector<size_t> prevBound(n_part, 0);
  std::vector<std::vector<double> > torques (n_part, std::vector<double>(3, default_value));
  std::vector<std::vector<double> > forces (n_part, std::vector<double>(3, default_value));
  std::vector<std::vector<double> > randomComponent (n_part, std::vector<double>(6, default_value));
  std::vector<std::vector<int> >   connectivity (n_part, std::vector<int>  (n_part, default_value));
  std::vector<int> stapleNumbers(n_part, default_value);
  std::vector<bool> isCrossover(n_part, false);
  std::vector<std::vector<int> > SM (n_part, std::vector<int> (n_part, default_value));
  std::vector<std::vector<int> >   staple_connections (n_part, std::vector<int>  (n_part, default_value));
  std::vector<position3D<double> > r(n_part);
  std::vector<int> belongsTo(n_part,  -1);
  std::vector<int> StrandNumber(n_part, -1);
  // CREATE COPIES OF FORCE, TORQUE, AND R FOR RUNGE KUTTA ALGORITHM
  std::vector<std::vector<double> > forces_forw = forces;
  std::vector<std::vector<double> > torques_forw = torques;
  std::vector<position3D<double> > r_proposed = r;
  std::vector<VerletList> neighbors;

  // RANDOM NUMBER GENERATORS
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::normal_distribution<double> distribution(0,1);
  srand(seed);
  size_t n_bound_staples = 0;
  load_file(&inFile[0], SM, staple_connections, stapleNumbers, connectivity, n_bonds, n_staple_seqs, isCrossover, StrandNumber, n_nucleotides_per_bead, CIRCULAR_SCAFFOLD);
  std::cout << "LOADED " << inFile << ".\n";

  if (CIRCULAR_SCAFFOLD && StrandNumber[0] != StrandNumber[n_scaf-1]){
    isCrossover[0]=true;
    isCrossover[n_scaf-1]=true;
  }
  
  for (size_t i = 0; i < n_part;i++){
    for (size_t j = 0; j < n_part; j++){
      if (SM[i][j]){
	belongsTo[i]=j;
	belongsTo[j]=i;
      }
    }
  }

  makeConnectivityMatrix(n_scaf, connectivity);
  for(size_t i = 0; i < n_part; i++){
    for (size_t j = 0; j < n_part; j++){
      if (SM[i][j]){
	belongsTo[i]=j;
      }
    }
  }  
  
  std::ofstream metadatafile;
  std::string metadatafilename;
  metadatafilename.append("simulation_metadata.dat");
  metadatafile.open(metadatafilename);  
  metadatafile << "Simulation Metadata\n";
  metadatafile << "Input JSON: " << inFile << "\n";
  metadatafile << "Output file: " << outFile << "\n";
  metadatafile << "Steps Per Frame: " << stepsPerFrame << "\n";
  metadatafile << "Simulation Length (steps): " << lsim << "\n";
  metadatafile << "dt: " << dt << " nanoseconds \n";
  metadatafile << "Printing to stdout every " << print_to_stdout_every << " steps\n";
  metadatafile << "Starting temp: " << temp << " K\n";
  metadatafile << "Final temp: " << final_temp << " K\n";
  metadatafile << "Using annealing rate of " << annealing_rate << " K / step\n";
  metadatafile << "Outputting last conf to last_conf file: " << lastconfname << "\n";
  metadatafile << "Circular scaffold is turned " << (CIRCULAR_SCAFFOLD ? "on" : "off") << "\n"; 
  metadatafile << "Forced binding is turned " << (FORCED_BINDING ? "on" : "off") << "\n";
  metadatafile << "Periodic boundaries are turned " << (pbc ? "on" : "off") << "\n";
  metadatafile << "Identified " << n_scaf << " scaffold beads and " << n_stap << " staple beads for a total of " << n_part << " beads\n"; 
  metadatafile << "Using hybridization strength of " << binding_energy_kcal_mol << " kcal/mol\n";
  metadatafile << "Using hybridization cutoff of " << binding_distance_cutoff << " nm\n";
  metadatafile.close();
  std::ofstream PerfectComplements;
  std::string PerfectComplementsName;
  PerfectComplementsName.append("complements.dat");
  PerfectComplements.open(PerfectComplementsName);
  PerfectComplements << "ScaffoldParticle PairedStapleParticle StrandNumber\n";
  for (size_t i = 0; i < n_scaf; i++){
    PerfectComplements << i << " " << belongsTo[i] << "\n";
  }
  PerfectComplements.close();
  std::string OUT_TYPE = "LAMMPS";
  ////////////////////////////////////////////////////////////////////
  // PRINT BASIC SIMULATION PARAMS //////////////////////////////////
  //////////////////////////////////////////////////////////////////
  std::cout << "Gamma_trans = " << gamma_trans <<".\n";
  std::cout << "Viscosity = " << visc << ".\n";
  std::cout << "Temp = " << temp << ".\n";
  std::cout << "Gamma_rot = " << gamma_rot << ".\n";
  std::cout << "Using seed = " << seed << ".\n";
  std::cout << "Initializing particles...\n";
  //////////////////////////////////////////////////////////////////
  // INIT /////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  if (!RESTART){
    simInit(r, n_scaf, simbox, beadDiameter, n_part, beadAxialSeparation, beadRadialSeparation, staple_connections, SM);
  } else {
    //load_frame(lastconfname,xboxsize,yboxsize,zboxsize,NATOM,STEPNUMBER,r);
  }
  if (RESTART){
    simbox.dimensions.x = xboxsize;
    simbox.dimensions.y = yboxsize;
    simbox.dimensions.z = zboxsize;
  }

  int n_bonds = 0;
  for (size_t i = 0; i < n_part; i++){
    for (size_t j = 0; j < n_part; j++){
      if (staple_connections[i][j]) { // TYPE 1: STAPLE BACKBONE
        n_bonds++;
      }
      if (connectivity[i][j]) { // TYPE 3: SCAFFOLD BACKBONE
        n_bonds++;
      }
    }
  }


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
  topf << n_bonds << " bonds\n";
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
    topf << "   " << i+1 << "   0   ";
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
      if (connectivity[i][j]) { // TYPE 1: SCAFFOLD BACKBONE
        bond_idx++;
        topf << "   " << bond_idx << "   1   " << i+1 << "   " << j+1 << "\n";
      }
      if (staple_connections[i][j]) { // TYPE 2: STAPLE BACKBONE
        bond_idx++;
        topf << "   " << bond_idx << "   2   " << i+1 << "   " << j+1 << "\n";
      }
    }
  }
  topf.close();
  std::cout << "printed topology file\n";
  /////////////////////////////////////////////////////////////////
  //// SIMULATION STARTS HERE ////////////////////////////////////
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
  std::cout << "Initialized clocks and reset simulation timer\n";
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
    ForceLoopsThreads.push_back(std::thread(calculatePairForces,std::ref(forces),std::ref(r),std::ref(simbox),std::ref(pbc),std::ref(isCrossover),std::ref(beadAxialSeparation),std::ref(beadRadialSeparation),std::ref(epsilon),std::ref(sigma),std::ref(n_part),std::ref(n_scaf),std::ref(n_stap),std::ref(belongsTo),std::ref(isBound),std::ref(staple_connections), std::ref(r_cut),std::ref(CIRCULAR_SCAFFOLD),std::ref(btout),std::ref(t),std::ref(n_bound_staples),true,std::ref(StrandNumber),std::ref(fbtOut),std::ref(prevBound),std::ref(binding_energy_kcal_mol),std::ref(binding_distance_cutoff),std::ref(FORCED_BINDING),std::ref(mutexes), std::ref(neighbors), NTHREADS));
    ForceLoopsThreads.push_back(std::thread(calculateBendingForces,n_scaf,std::ref(isBound),std::ref(r), std::ref(forces), dsdna_lp, k_B*temp, beadAxialSeparation, std::ref(stapleNumbers), std::ref(isCrossover), std::ref(simbox), pbc));
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
    RungeKuttaThreads.push_back(std::thread(calculatePairForces,std::ref(forces),std::ref(r_proposed),std::ref(simbox),std::ref(pbc),std::ref(isCrossover),std::ref(beadAxialSeparation),std::ref(beadRadialSeparation),std::ref(epsilon),std::ref(sigma),std::ref(n_part),std::ref(n_scaf),std::ref(n_stap),std::ref(belongsTo),std::ref(isBound),std::ref(staple_connections), std::ref(r_cut),std::ref(CIRCULAR_SCAFFOLD),std::ref(btout),std::ref(t),std::ref(n_bound_staples),true,std::ref(StrandNumber),std::ref(fbtOut),std::ref(prevBound),std::ref(binding_energy_kcal_mol),std::ref(binding_distance_cutoff),std::ref(FORCED_BINDING),std::ref(mutexes), std::ref(neighbors), NTHREADS));
    RungeKuttaThreads.push_back(std::thread(calculateBendingForces,n_scaf,std::ref(isBound),std::ref(r_proposed), std::ref(forces), dsdna_lp, k_B*temp, beadAxialSeparation, std::ref(stapleNumbers), std::ref(isCrossover), std::ref(simbox), pbc));
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
    //Check that forces are balanced.
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
    // Real integration with averaged forces
    integrateMotion(n_part, forces, torques,  dt, r, r, randomComponent, verbose, pbc, simbox, n_scaf, gamma_trans, t, isCrossover, isBound);
    stop_frame = std::chrono::high_resolution_clock::now();
    duration_frame = std::chrono::duration_cast<std::chrono::microseconds>(stop_frame-start_frame);
    if (t%1001==0){
      std::cout << "Complete Timestep for step number" << t << " = " << duration_frame.count() << " microseconds.\n";
      std::cout << "Current temperature: " << temp << " Kelvin.\n";
      std::cout << "Currently " << n_bound_staples << " staples are bound\n";
    }
    // Adjust temp and fluid properties
    if (temp > final_temp){
      std::cout << "t = " << t*dt*1E-9*5200 << " / " << lsim*dt*1E-9*5200 << ".\n";
      temp -= annealing_rate;
      tm273 = temp - 273.15;
      visc = (3.245157366681122E-11 * pow(tm273,4) - 9.061289916246450E-09 * pow(tm273,3) + 9.845093457119180E-07 * pow(tm273,2) - 5.521101038709857E-05 * tm273 + 1.778370453327189E-03) * 1E3;
      gamma_trans = 4 * _PI * visc * (beadAxialSeparation+beadRadialSeparation)/4;
      Er_linker = 4 * _PI * visc * pow(beadAxialSeparation/2,2) * beadAxialSeparation;
    }
    t++;
  } // end of time loop
  trajectory.close();
  //torqueOutput.close();
  //forceOutput.close();
  btout.close();
  return EXIT_SUCCESS;
}
//
