#include "headers.hpp"
#include <cassert>
#include <iostream>
#include <ostream>
#include <istream>
#include <random>
#include <stdlib.h>
#include <stdio.h>

#ifndef PARAMETERS
#define PARAMETERS

#define _PI 3.14159265358979
#define DEIGEN_STACK_ALLOCATION_LIMIT 0

//SIMULATION UNITS FOR THIS PROGRAM:
// DISTANCE: nm
// TIME: ns
// FORCE: pN
// TEMPERATURE: K

// Parameters. Include in all source files.
bool verbose = false;
bool preBranch = false;
bool pbc = false;
bool timed = false;

bool output_TCL = false;
bool output_OVITO = false;
bool output_LAMMPS = true;
bool CIRCULAR_SCAFFOLD = false;
bool delay_binding = false;





// SQUARE OR HEX 
double n_nucleotides_per_bead = 8;

// BINDING PARAMETERS
double binding_energy_kcal_mol = 8; //kcal per mol
double binding_distance_cutoff = 2; //nm



//432
//402 for old 6hb, 410 for rerouted
const size_t n_scaf = 48;//224;//432;//224;//64;//432;// 48;//48;//28;// 48; // 512/8;
const size_t n_stap = 48;//224-8+1;//432;//224-8+1;//0;//432;//48; //28;//48; //512/8;

// Qinyi 1008 scaffold 972 staples
//const size_t n_scaf = 1008;
//const size_t n_stap = 972;
const size_t n_part = n_scaf + n_stap;
size_t n_staple_seqs = 0;
size_t n_bonds = 0;

// Shrinking periodic box condition
double shrink_rate = 0.00; //nm/step
double simbox_final_size_ratio = 1; //unitless

// Time-based parameters
const size_t lsim = 5E8; //steps
const double dt = .008; //nanoseconds
const size_t stepsPerFrame = 10000;
const size_t print_to_stdout_every = 20000;

// Environment  
const double temp = 300.0; //Kelvin
double tm273 = temp - 273.15E0;
double visc = 1.77837;//  //(3.245157366681122E-11 * pow(tm273,4) - 9.061289916246450E-09 * pow(tm273,3) + 9.845093457119180E-07 * pow(tm273,2) - 5.521101038709857E-05 * tm273 + 1.778370453327189E-03) * 1E3;;

//std::cout << "running simulation for temp = " << temp << "\n" << "visc = " << visc << "\n";

// Particle interactions
const double beadDiameter = 3; //nm
const double beadAxialSeparation = 3.15; // 2.656;// * 7/8;// * 7/8; //nm //MODIFIED, MAKE SURE THIS IS RIGHT
const double beadRadialSeparation = 2.8; //nm

const double sigma = pow(beadRadialSeparation,0.833);
const double epsilon = 6.96; //pN*nm


const double r_cut = 15; //nm
const double l_k = .004; //inverse spring constant of bead-chain
const double crossover_stiffness = 30; //N/m
const double dsdna_lp = 50; //nm


// THIS IS WHERE WE PRINT LAMMPS STYLE TOPOLOGY FILE
std::string topology_filename = "DirectSimulationTopology";

// Simulation box, same size as 6hb
double CubicBoxSize = ((2.656 * n_scaf) + 2 * r_cut) / 4  ;
simBox<double> simbox(CubicBoxSize, CubicBoxSize, CubicBoxSize);

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

//std::vector<double> Alpha (n_part-1, default_size_t);
//std::vector<double> Beta (n_part-1, default_size_t);
//std::vector<double> Gamma (n_part-1, default_size_t);

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

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
std::normal_distribution<double> distribution(0,1);
size_t n_bound_staples = 0;

#endif
