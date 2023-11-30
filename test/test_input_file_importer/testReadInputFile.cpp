#include "../../src/readInputFile.hpp"
using namespace std;

int main (void){
  vector<pair<string,string> > inputParams;
  std::string inputFileName = "testInput.txt";
  int success = readInputFile(inputFileName,inputParams);
  for (size_t i = 0; i < inputParams.size(); ++i){
    std::cout << inputParams[i].first << "=" << inputParams[i].second << "\n";
  }
  return EXIT_SUCCESS;
}
  
