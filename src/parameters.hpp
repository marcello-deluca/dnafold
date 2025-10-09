#include "headers.hpp"
#include <cassert>
#include <iostream>
#include <fstream>
#include <ostream>
#include <istream>
#include <random>
#include <stdlib.h>
#include <stdio.h>
#include <thread>
#include <mutex>
//#include "loadFrame.hpp"
#include "VerletList.hpp"
#include "MakeVerletList.hpp"
#ifndef PARAMETERS
#define PARAMETERS
// SIMULATION UNITS FOR THIS PROGRAM:
// DISTANCE: nm
// TIME: ns
// FORCE: pN
// TEMPERATURE: K
// ENERGY: pN-nm
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
double n_nucleotides_per_bead = 8;
// BINDING PARAMETERS
double binding_energy_kcal_mol = 10; //kcal per mol
double binding_distance_cutoff = 2; //nm
bool FORCED_BINDING = true;
// NUMBER OF SCAFFOLD AND STAPLE IN DESIGN
size_t n_scaf = 0;
size_t n_stap = 0;
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
const double beadAxialSeparation = 2.725; //nm
const double beadRadialSeparation = 2.4; //nm
const double sigma = beadRadialSeparation/1.12;
const double epsilon = 6.96; //pN*nm
const double r_cut = 15; //nm
const double l_k = .004; //inverse spring constant of bead-chain
const double crossover_stiffness = 30; //N/m
const double dsdna_lp = 50; //nm
// SIMULATION BOX SIZE
double CubicBoxSize = 200;
// SHRINKING BOX CONDITION
double shrink_rate = 1E-5; //nm/step
double simbox_final_size_ratio = .5; //unitless
size_t n_staple_seqs = 0;
size_t n_bonds = 0;
// Initializing simluation variables double sigma;
double dr;
double k_B = .0138; //pN*nm / K
double gamma_trans = 6 * _PI * visc * (beadAxialSeparation+beadRadialSeparation)/4; // 4πηR   gamma units: pN / nm * ns
double gamma_rot  = 8 * _PI * visc * pow(beadDiameter/2,3); // 8πηR^3
double Er_linker = 4 * _PI * visc * pow(beadAxialSeparation/2,2) * beadAxialSeparation;
double r_var = 2 * dt * k_B * temp / gamma_rot;
double stretch;
double e_dist;
double default_value = 0;
size_t default_size_t = 0;
#endif
