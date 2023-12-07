#ifndef PARSE_JSON_H
#define PARSE_JSON_H

#include<string>
#include"nlohmann/json.hpp"
#include <vector>
#include <iostream>

void parse_json (char * ifname, size_t & n_scaf, size_t & n_stap, std::vector<std::vector<int> > & Scaffolds, std::vector<std::vector<int> > & Staples){
  FILE * f = fopen(ifname, "r");
  if (f==NULL){
    std::cerr << "Invalid file input\n";
  }
  std::string jsonstring;
  char c;
  while ((c=fgetc(f))!=EOF){
    jsonstring += c;
  }
  auto j = nlohmann::json::parse(jsonstring);

  std::vector<int> ThreePrimeScaffoldEnd;
  std::vector<std::vector<int> > ThreePrimeStapleEnds;
  int curr_strand;
  bool verbose = false;
  int prev_strand=-1;
  bool first_item = true;
  std::cout << "Parsing JSON file..." << std::endl;
  // Part 1: parse through the JSON file, and make lists called "Scaffolds" and "staples" containing the scaffolds and staples, their strand and index numbers (added custom), and their 3' and 5' neighbors 
  // formatted: Staples: vector<vector<int> >, e.g. Staples[i]={3,27, 3, 26, 3, 28 } would be a staple on strand 3 ind 27 with neighbors at 26 and 28
  for (auto& el : j["vstrands"].items()){
    if (verbose){
      std::cout << el.key() << " : " << el.value() << '\n' << '\n';
    }
    for (auto& el2 : el.value().items()) {
      if (el2.key() == "num"){
	      curr_strand = el2.value();
      }
      if (el2.key() == "scaf"){
      if (verbose){
        std::cout << "scaffolds:" << '\n';
      }
      for (auto& el3 : el2.value().items()){
        std::vector<int> current_scaffold;
        
        current_scaffold.push_back(curr_strand);
        current_scaffold.push_back(std::stoi(el3.key()));
        for (auto& el4 : el3.value().items()){
          if (verbose){
            std::cout << el4.value();
          }
          std::string strval = el4.value().dump();
          current_scaffold.push_back(std::stoi(strval));
        }
        for (size_t pq = 0; pq < current_scaffold.size(); pq++){
          if (verbose){  
            std::cout << current_scaffold[pq];
          }
          if (pq < current_scaffold.size()-1){
            if (verbose){
        std::cout << " , ";
            }
        }
        }
        if (current_scaffold[2] == -1 && current_scaffold [3] == -1 && current_scaffold[4]!=-1 &&current_scaffold[5]!=-1){
          ThreePrimeScaffoldEnd = current_scaffold; 
        }
        Scaffolds.push_back(current_scaffold);
        if (verbose){
          std::cout << "   " << el3 << '\n';
        }
      }
    }
      else if (el2.key() == "stap"){
        for (auto& el3 : el2.value().items()){
          std::vector<int> current_staple;
          current_staple.push_back(curr_strand);
          current_staple.push_back(std::stoi(el3.key()));
          for (auto& el4 : el3.value().items()){
            std::string strval = el4.value().dump();
            current_staple.push_back(std::stoi(strval));
          }
          for (size_t pq = 0; pq < current_staple.size();pq++){
            if (pq < current_staple.size()-1){
            }
          }
          Staples.push_back(current_staple);
          if (current_staple[2]==-1 && current_staple[3]==-1 && current_staple[4]!=-1)
          ThreePrimeStapleEnds.push_back(current_staple);
        }
      }
      else {
      }
    }
    prev_strand = curr_strand;
    first_item=false;
  }
  size_t n_staple_nucleotides = 0;
  size_t n_scaffold_nucleotides = 0;
  for (size_t i = 0; i < Scaffolds.size(); ++i){
    if (Scaffolds[i][2]!=-1 || Scaffolds[i][4]!=-1){
      n_scaffold_nucleotides++;
    }
  }
  for (size_t i = 0; i < Staples.size(); ++i){
    if (Staples[i][2]!=-1 || Staples[i][4]!=-1){
      n_staple_nucleotides++;
    }
  }
  std::cout << "found " << n_scaffold_nucleotides << " scaffold nucleotides and " << n_staple_nucleotides << " staple nucleotides\n";
  n_stap = n_staple_nucleotides/8;
  n_scaf = n_scaffold_nucleotides/8;
  std::cout << "Finished parsing JSON file.\n";
}
#endif
