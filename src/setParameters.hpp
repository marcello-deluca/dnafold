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
void setParameters (size_t & stepsPerFrame, size_t & lsim, double & dt, size_t & print_to_stdout_every, size_t & stepsPerNeighborListRefresh, double & temp, double & final_temp, double & annealing_rate, std::string & lastconfname, bool & CIRCULAR_SCAFFOLD, bool & FORCED_BINDING, size_t NTHREADS, int & reportTimings, string & inFile, string & outFile, vector<pair<string,string> > & inputs ){
  for (size_t i = 0; i < inputs.size(); ++i){
    std::string val = inputs[i].second;
    switch (inputs[i]){
    case "stepsPerFrame":
      stepsPerFrame = (size_t) stoi(val);
      break;
    case "lsim":
      lsim = (size_t) stoi(val);
      break;
    case "dt":
      dt = stod(val);
      break;
    case "print_to_stdout_every":
      print_to_stdout_every = (size_t) stoi(val);
      break;
    case "stepsPerNeighborListRefresh":
      stepsPerNeighborListRefresh = (size_t) stoi(val);
      break;
    case "temp":
      temp = stod(val);
      break;
    case "final_temp":
      final_temp = stod(val);
      break;
    case "annealing_rate":
      annealing_rate = stod(val);
      break;
    case "lastconfname":
      lastconfname = val;
      break;
    case "CIRCULAR_SCAFFOLD":
      CIRCULAR_SCAFFOLD = evaluateTrueFalseStatement(val);
      break;
    case "FORCED_BINDING":
      FORCED_BINDING = evaluateTrueFalseStatement(val);
      break;
    case "NTHREADS":
      NTHREADS = (size_t) stoi(val);
      break;
    case "ReportTimings":
      ReportTimings = stoi(val);
      break;
    case "inFile":
      inFile = val;
      break;
    case "outFile":
      outFile = val;q
      break;
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
