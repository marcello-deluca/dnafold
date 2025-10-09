#ifndef SETPARAMETERS_H
#define SETPARAMETERS_H
#include <vector>
using namespace std;
bool evaluateTrueFalseStatement(std::string in){
  if (in=="true"|| in=="1"){
    return true;
  } else if (in=="false" || in=="0"){
    return false;
  } else {
    std::cerr << "bad true/false output\n";
    exit(0);
  }
}
void setParameters (size_t & stepsPerFrame, size_t & lsim, double & dt, size_t & print_to_stdout_every, size_t & stepsPerNeighborListRefresh, double & temp, double & final_temp, double & annealing_rate, std::string & lastconfname, bool & CIRCULAR_SCAFFOLD, bool & FORCED_BINDING, size_t & NTHREADS, int & reportTimings, bool & pbc, string & inFile, string & outFile, double & binding_energy_kcal_mol, double & CubicBoxSize, double & shrink_rate, double & simbox_final_size_ratio, vector<pair<string,string> > & inputs){
  for (size_t i = 0; i < inputs.size(); ++i){
    std::string key = inputs[i].first;
    std::string val = inputs[i].second;
    if (key ==  "stepsPerFrame"){
      stepsPerFrame = (size_t) stoi(val);
      std::cout << "Printing to trajectory file every " << stepsPerFrame << " frames\n";
    }
    if (key == "lsim"){
      lsim = (size_t) stod(val);
      std::cout << "simulation length: " << lsim << " steps\n";
    }
    if (key == "dt"){
      dt = stod(val);
      std::cout << "using timestep of " << dt << " nanoseconds\n";
    }
    if (key == "print_to_stdout_every"){
      print_to_stdout_every = (size_t) stoi(val);
      std::cout << "Printing to standard output every " << print_to_stdout_every << " steps\n";
    }
    if (key == "stepsPerNeighborListRefresh"){
      stepsPerNeighborListRefresh = (size_t) stoi(val);
      std::cout << "Refreshing neighbor list every " << stepsPerNeighborListRefresh << " steps\n";
    }
    if (key == "temp"){
      temp = stod(val);
      std::cout << "Starting temp: " << temp << " K\n";
    }
    if (key == "final_temp"){
      final_temp = stod(val);
      std::cout << "Final temp: " << final_temp << " K\n";
    }
    if (key == "annealing_rate"){
      annealing_rate = stod(val);
      std::cout << "Using annealing rate of " << annealing_rate << " K / step\n";
    }
    if (key == "lastconfname"){
      lastconfname = val;
      std::cout << "Outputting last conf to last_conf file: " << lastconfname << "\n";
    }
    if (key == "CIRCULAR_SCAFFOLD"){
      CIRCULAR_SCAFFOLD = evaluateTrueFalseStatement(val);
    }
    if (key == "FORCED_BINDING"){
      FORCED_BINDING = evaluateTrueFalseStatement(val);
    }
    if (key == "NTHREADS"){
      NTHREADS = (size_t) stoi(val);
      std::cout << "Carrying out force calculations using " << NTHREADS << " threads\n";
    }
    if (key == "reportTimings"){
      reportTimings = stoi(val);
    }
    if (key == "inFile"){
      inFile = val;
      std::cout << "Using " << inFile << " as reference design file\n";
    }
    if (key == "outFile"){
      outFile = val;
      std::cout << "Outputting trajectory to " << outFile << "\n";
    }
    if (key == "pbc"){
      pbc = evaluateTrueFalseStatement(val);
      if (pbc){
	std::cout << "using periodic boundary conditions\n";
      } else {
	std::cout << "not using periodic boundary conditions\n";
      }
    }
    if (key == "BoxSize"){
      CubicBoxSize = stod(val);
      std::cout << "Using starting box size of " << CubicBoxSize << "\n";
    }
    if (key == "ShrinkRate"){
      shrink_rate = stod(val);
      std::cout << "Shrinking box at rate of " << shrink_rate << " nm/frame\n";
    }
    if (key == "FinalSizeRatio"){
      simbox_final_size_ratio = stod(val);
      std::cout << "Final size ratio set to " << simbox_final_size_ratio << " for a final box size of " << simbox_final_size_ratio * shrink_rate << " nm\n";
    }
    if (key == "bindStrength"){
      binding_energy_kcal_mol = stod(val);
      std::cout << "using bead binding enthalpy of " << binding_energy_kcal_mol << " kcal/mol\n";
    }
  }
}

/*
   parameters that we can adjust in the input file:
   stepsPerFrame (size_t)
   lsim (size_t)
   dt (double)
   print_to_stdout_every (size_t)
   stepsPerNeighborListRefresh (size_t)
   temp (double)
   final_temp (double)
   annealing_rate (double)
   fluid (not currently implemented)
   lastconfname (std::string)
   CIRCULAR_SCAFFOLD (bool)
   FORCED_BINDING (bool)
   NTHREADS (size_t)
   ReportTimings (unimplemented, 0 for false or size_t for every ReportTimings steps)
   inFile (string)
   outFile (string)
*/

#endif
