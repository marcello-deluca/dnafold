#ifndef READINPUTFILE_H
#define READINPUTFILE_H

#include <iostream>
#include <fstream>
#include <vector>

/*
  Pseudocode:
  - Input should be structured "input_variable=value" or "input_variable = value"  or any combination of these / any number of spaces.
  - Parse the input into a vector of tuples
     • Grab each string for the line of text
     • Remove special characters like newline
     • Parse out the first statement (variable in question) and second statement (value to set that  variable to)
     • Possible bad input is inputting a double for an int or size_t-typed variable. Need to find a way to deal with this.
  - Go through a switch/case for each input to set its value
  - If the user inputs a case that is not in the list of options, terminate the program and let them know that they made a mistake.
  - Once everything is done and all variables are set, output some text in stdout to let the user know that each variable has been set to the requested value in list form for troubleshooting
  - Comments and weird stuff not considered
*/

void cleanString(std::string & in){
  bool comment = false;
  size_t cLoc = in.find("#");
  if (cLoc!=std::string::npos){
    in=in.substr(0,cLoc-1);
  }
  size_t spLocation = in.find(" ");
  while (spLocation!=std::string::npos){
    in.erase(spLocation,1);
    spLocation=in.find(" ");
  }
}

void cleanInput(std::string & LHS, std::string & RHS){
  cleanString(LHS);
  cleanString(RHS);
}

std::pair<std::string, std::string> parseLine(std::string inputLine){
  size_t eqLocation = inputLine.find("=");
  std::string LHS, RHS;
  if (eqLocation != std::string::npos) {
    LHS = inputLine.substr(0, eqLocation);
    RHS = inputLine.substr(eqLocation+1,inputLine.size());
    cleanInput(LHS,RHS);
  } else {
    std::cerr << "bad input\n";
    exit(0);
  }
  std::pair<std::string,std::string> out(LHS,RHS);
  return out;
}
  
int readInputFile(std::string fileLocation, std::vector<std::pair<std::string, std::string> > & InputParams){
  std::string line;
  std::ifstream inFile(fileLocation);
  if (inFile.is_open()){
      while (getline(inFile,line)){
	std::cout << line <<"\n";
	if(line.find("#")==0 || line==""){} else {
	  InputParams.push_back(parseLine(line));
	}
      }
      for (size_t i = 0; i < InputParams.size(); ++i){
	std::cout << "Parameter " << i << " values: (first) " << InputParams[i].first << "; (second) " << InputParams[i].second << "\n";
      }
  } else {
    std::cerr << "bad input\n";
    return -1;
  }
  inFile.close();
  return 1;
}

#endif
