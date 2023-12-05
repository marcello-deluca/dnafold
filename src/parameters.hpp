#include "headers.hpp"
#include <cassert>
#include <iostream>
#include <fstream>
#include <ostream>
#include <istream>
#include <random>
#include <stdlib.h>
#include <stdio.h>
//#include "omp.h"
#include <execution>
#include <thread>
#include <mutex>
#include"loadFrame.hpp"
#include "VerletList.hpp"
#include "MakeVerletList.hpp"

#ifndef PARAMETERS
#define PARAMETERS
#define _PI 3.14159265358979
#define DEIGEN_STACK_ALLOCATION_LIMIT 0

// SIMULATION UNITS FOR THIS PROGRAM:
// DISTANCE: nm
// TIME: ns
// FORCE: pN
// TEMPERATURE: K
// ENERGY: pN-nm

// Parameters. Include in all source files.
bool verbose = false;
bool preBranch = false;
bool timed = false;
bool delay_binding = false;
bool output_TCL = false;
bool output_OVITO = false;
bool output_LAMMPS = true;

// PERIODIC BOUNDARY CONDITIONS
bool pbc = false;

// IS SCAFFOLD CIRCULAR?
bool CIRCULAR_SCAFFOLD = false;

// SQUARE OR HEX 
double n_nucleotides_per_bead = 8; //8 = square, 7 = hex

// BINDING PARAMETERS
double binding_energy_kcal_mol = 10; //kcal per mol
double binding_distance_cutoff = 2; //nm
bool FORCED_BINDING = true;

// NUMBER OF SCAFFOLD AND STAPLE IN DESIGN
const size_t n_scaf = 224;// 48;//432;//224;//64;//432;// 48;//48;//28;// 48; // 512/8;
const size_t n_stap = 224;// 48;//432;//224-8+1;//0;//432;//48; //28;//48; //512/8;
std::vector<std::mutex> mutexes(n_scaf+n_stap);

// TIME BASED PARAMETERS
size_t lsim = 1E9; //steps
double dt = .01; //nanoseconds
size_t stepsPerFrame = 10000;
size_t print_to_stdout_every = 20000;
size_t stepsPerNeighborListRefresh = 1000;


// ENVIRONMENT
double temp = 300.0; //Kelvin
double tm273 = temp - 273.15E0;
double visc = (3.245157366681122E-11 * pow(tm273,4) - 9.061289916246450E-09 * pow(tm273,3) + 9.845093457119180E-07 * pow(tm273,2) - 5.521101038709857E-05 * tm273 + 1.778370453327189E-03) * 1E3;

// ANNEALING
double final_temp = 300.0;
double annealing_rate = 0.1E-5; //Annealing rate in K/step

// PARTICLE INTERACTIONS
const double beadDiameter = 2.7; //nm
const double beadAxialSeparation = 2.725; // 2.656;// * 7/8;// * 7/8; //nm //MODIFIED, MAKE SURE THIS IS RIGHT
const double beadRadialSeparation = 2.4; //nm
const double sigma = pow(beadRadialSeparation,0.833);
const double epsilon = 6.96; //pN*nm
const double r_cut = 15; //nm
const double l_k = .004; //inverse spring constant of bead-chain
const double crossover_stiffness = 30; //N/m
const double dsdna_lp = 50; //nm

// SIMULATION BOX SIZE
double CubicBoxSize = 200;
simBox<double> simbox(CubicBoxSize, CubicBoxSize, CubicBoxSize);

// SHRINKING BOX CONDITION
double shrink_rate = 1E-5; //nm/step
double simbox_final_size_ratio = .5; //unitless

const size_t n_part = n_scaf + n_stap;
size_t n_staple_seqs = 0;
size_t n_bonds = 0;

// THIS IS WHERE WE PRINT LAMMPS STYLE TOPOLOGY FILE
//std::string topology_filename = "DirectSimulationTopology";



// Initializing simluation variables double sigma;
double dr;
double k_B = .0138; //       pN*nm / K

//modified to go faster
double gamma_trans = 6 * _PI * visc * (beadAxialSeparation+beadRadialSeparation)/4; // 4πηR   gamma units: pN / nm * ns
double gamma_rot  = 8 * _PI * visc * pow(beadDiameter/2,3); // 8πηR^3
double Er_linker = 4 * _PI * visc * pow(beadAxialSeparation/2,2) * beadAxialSeparation;

// CHANGED TO gamma_rot from er_linker 
double r_var = 2 * dt * k_B * temp / gamma_rot;
double stretch;
double e_dist;
double default_value = 0;
size_t default_size_t = 0;

// INITIALIZE VECTORS
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
size_t n_bound_staples = 0;

// RESTART SPECS AND VARIABLES
bool RESTART = false;
std::string lastconfname = "last_conf.dat";
double xboxsize, yboxsize, zboxsize;
int NATOM;
size_t STEPNUMBER = 0;


#endif
