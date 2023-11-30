//-std=c++11A
#include "potentials.hpp"
#include "forceUtils.hpp"
#include "headers.hpp"
#include "simBox.hpp"
#include "position3D.hpp"
#include "Interaction.hpp"
#include <thread>
#include <mutex>
#include <cmath>
std::mutex mtx;
void calculateIndividualParticleForces(size_t firstIndex, size_t secondIndex, std::vector<std::vector<double> > & forces, std::vector<position3D<double> > & r, simBox<double> & sb, bool pbc, std::vector<bool> & isCrossover, double beadAxialSeparation, double beadRadialSeparation, double epsilon, double sigma, const size_t n_part, const size_t n_scaf, const size_t n_stap, std::vector<int> & belongsTo, std::vector<size_t> & isBound, std::vector<std::vector<int> > & staple_connections, double RCUT,  bool CIRCULAR_SCAFFOLD, std::ofstream & btout, size_t t, size_t & n_bound_staples, bool PRE_RUNGE_KUTTA, std::vector<int> & StrandNumber, std::ofstream & fbtOut, std::vector<size_t> & prevBound, double bind_energy_kcal_mol, double BIND_CUTOFF, bool FORCED_BINDING, std::vector<std::mutex> & mtx, std::vector<Interaction> & InteractionList){
  bool verbose = false;
  bool excVerbose = false;
  double k_stretch = 152;
  double RXI, RYI, RZI;
  double RXIJ = 0;
  double RYIJ = 0;
  double RZIJ = 0;
  double DIST = 0;
  double FORCE;
  int IDX;
  double ssdna_dist = beadAxialSeparation;
  double ssdna_k = k_stretch;
  double eqdist = 0;
  double bind_energy = bind_energy_kcal_mol * 6.97; //convert to pN*nm
  double bind_dist = BIND_CUTOFF; //nanometers
  double bind_force = bind_energy/bind_dist; //piconewton
  double FORCED_BINDING_F = 10;
  double exc_vol_const =  6.97; // 6.97 pN*nm = 1kcal/mol
  size_t nScm1 = n_scaf-1;
  double dx, dy, dz;

  for (size_t i = firstIndex; i <= secondIndex; ++i){
    Interaction vli = InteractionList[i];
    int FirstParticleIndex = vli.First;
    int SecondParticleIndex = vli.Second;
    RXI = r[FirstParticleIndex].x;
    RYI = r[FirstParticleIndex].y;
    RZI = r[FirstParticleIndex].z;
    size_t sz = vli.InteractionList.size();
    //STAPLE-SCAFFOLD BINDING POTENTIAL
    if (FirstParticleIndex < n_scaf){
      
    } else { //SECOND HALF OF MAIN LOOP: STAPLES
    } //close "else"


    // EXCLUDED VOLUME INTERACTIONS
    for (size_t nidx = 0; nidx < sz; ++nidx){
      size_t idx = vli.NeighborList[nidx];
      if (i < n_scaf && idx >= n_scaf && belongsTo[i]!=idx){ // SCAFFOLD-STAPLE EXCLUDED VOLUME
	distances(pbc,sb,RXI,RYI,RZI,RXIJ,RYIJ,RZIJ,DIST,r,idx);
	if (DIST > RCUT) {
	} else {
	  FORCE = excvolwca(sigma, exc_vol_const/10,DIST);
	  //addPairForces(i,idx,forces,FORCE,RXIJ,RYIJ,RZIJ,DIST, mtx);
	}
      } else if (i < n_scaf && idx < n_scaf && abs((int)(i-idx))>1 ){ // SCAFFOLD SCAFFOLD EXCLUDED VOLUME
	distances(pbc,sb,RXI,RYI,RZI,RXIJ,RYIJ,RZIJ,DIST,r,idx);
	if (DIST > RCUT) {
	} else {
	  FORCE = excvolwca(sigma, exc_vol_const,DIST);
	  addPairForces(i,idx,forces,FORCE,RXIJ,RYIJ,RZIJ,DIST, mtx);
	}
      } else if (i > n_scaf && idx > n_scaf && abs((int)(i-idx))>1){
	distances(pbc,sb,RXI,RYI,RZI,RXIJ,RYIJ,RZIJ,DIST,r,idx);
	if (DIST > RCUT) {
	} else {
	  FORCE = excvolwca(sigma, exc_vol_const,DIST);
	  addPairForces(i,idx,forces,FORCE,RXIJ,RYIJ,RZIJ,DIST, mtx);
	}
      }
    }
  } // end of i iterator
}

void calculatePairForces(std::vector<std::vector<double> > & forces, std::vector<position3D<double> > & r, simBox<double> & sb, bool pbc, std::vector<bool> & isCrossover, double beadAxialSeparation, double beadRadialSeparation, double epsilon, double sigma, const size_t n_part, const size_t n_scaf, const size_t n_stap, std::vector<int> & belongsTo, std::vector<size_t> & isBound, std::vector<std::vector<int> > & staple_connections, double RCUT,  bool CIRCULAR_SCAFFOLD, std::ofstream & btout, size_t t, size_t & n_bound_staples, bool PRE_RUNGE_KUTTA, std::vector<int> & StrandNumber, std::ofstream & fbtOut, std::vector<size_t> & prevBound, double bind_energy_kcal_mol, double BIND_CUTOFF, bool FORCED_BINDING, std::vector<std::mutex> & mtx, std::vector<Interaction> & InteractionList){
    bool verbose = false;
    bool excVerbose = false;
    double k_stretch = 152;
    double RXI, RYI, RZI;
    double RXIJ = 0;
    double RYIJ = 0;
    double RZIJ = 0;
    double DIST = 0;
    double FORCE;
    int IDX;
    double ssdna_dist = beadAxialSeparation;
    double ssdna_k = k_stretch;
    double eqdist = 0;
    double bind_energy = bind_energy_kcal_mol * 6.97; //convert to pN*nm
    double bind_dist = BIND_CUTOFF; //nanometers
    double bind_force = bind_energy/bind_dist; //piconewton
    double FORCED_BINDING_F = 1;
    double exc_vol_const =  6.97; //pN*nm; 6.97 pN*nm = 1kcal/mol 
    size_t nScm1 = n_scaf-1;
    size_t NTHREADS = 2;

    // FIRST SECTION: ENFORCE CIRCULAR SCAFFOLD
    if (CIRCULAR_SCAFFOLD){
      if (StrandNumber[0]==StrandNumber[nScm1]){
	eqdist = beadAxialSeparation;
      } else {
	eqdist = beadRadialSeparation;
      }
      RXI = r[0].x;
      RYI = r[0].y;
      RZI = r[0].z;
      distances(pbc,sb,RXI,RYI,RZI,RXIJ,RYIJ,RZIJ,DIST,r,n_scaf-1);
      FORCE = simpleHarmonic(eqdist,DIST,k_stretch);
      addPairForces(0,n_scaf-1,forces,FORCE,RXIJ,RYIJ,RZIJ,DIST, mtx);
    }
    // MAIN LOOP: ITERATE THROUGH ALL PARTICLES
    n_interactions = InteractionList.size();
    size_t n_intx_per_thread = ceil(n_interactions / NTHREADS);
    std::vector<std::thread> threads;
    for (size_t i = 0; i < NTHREADS; ++i){
      size_t arg1 = i*n_intx_per_thread;
      size_t arg2 = (i+1)*n_intx_per_thread - 1;
      if (arg2 > n_part-1){
	std::cout << "weird duplicate force calculation problem\n";
	arg2 = n_part-1;
      }
      threads.push_back(std::thread(calculateIndividualParticleForces,arg1,arg2, std::ref(forces), std::ref(r), std::ref(sb), pbc, std::ref(isCrossover), beadAxialSeparation, beadRadialSeparation, epsilon, sigma, n_part, n_scaf, n_stap, std::ref(belongsTo), std::ref(isBound), std::ref(staple_connections), RCUT, CIRCULAR_SCAFFOLD, std::ref(btout), t, std::ref(n_bound_staples), PRE_RUNGE_KUTTA, std::ref(StrandNumber), std::ref(fbtOut), std::ref(prevBound), bind_energy_kcal_mol, BIND_CUTOFF, FORCED_BINDING,std::ref( mtx), std::ref(InteractionList)) );
    
    } //end "i" loop
    for( auto & thread : threads){
      while (!thread.joinable()){
      }
      thread.join();
    }
  } // Close function
