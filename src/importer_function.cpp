#include <cstdio>
#include <string>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include "nlohmann/json.hpp"
#include <vector>

using json = nlohmann::json;

//double n_nucleotides_per_bead = 8;



//////////////
//// FIND ////
//////////////

// This function takes a strand-index location and searches a list of beads for that location,
// then returns the index of that bead. 
int Find(int strand, int index, std::vector<std::vector<int> > & list){
  for (size_t i = 0; i < list.size();i++){
    if (list[i][0] == strand && list[i][1] == index){
      if (list[i][2]==-1 && list[i][3]==-1 && list[i][4]==-1 && list[i][5]==-1){
	//std::cout << "no scaffold at index " << strand << ", " << index << ".\n";
	return -1;
      }
      return i;
    }
  }
  std::cout << "//weird undefined case encountered" << std::endl;
  return -1;
}



///////////////////
//// LOAD FILE ////
///////////////////

// This program takes a single argument (a json file describing a cadnano design) and converts
// its contents into a format that can be simulated by the dna_bd model.
void load_file (char * ifname, std::vector<std::vector<int> > & SM_, std::vector<std::vector<int> > & staple_connections_, std::vector<int> & stapleNumbers_, std::vector<std::vector<int> > & connectivity_, size_t & n_bonds_, size_t & n_staple_seqs_, std::vector<bool> & isCrossover, std::vector<int> & StrandNumber, double n_nucleotides_per_bead, bool circular_scaffold){
  FILE * f = fopen(ifname, "r");
  if (f==NULL){
    std::cerr << "Invalid file input\n";
  }
  std::string jsonstring;
  char c;
  while ((c=fgetc(f))!=EOF){
    jsonstring += c;
  }
  auto j = json::parse(jsonstring);

  std::vector<std::vector<int> > Scaffolds;
  std::vector<std::vector<int> > Staples;
  std::vector<int> ThreePrimeScaffoldEnd;
  std::vector<std::vector<int> > ThreePrimeStapleEnds;
  int curr_strand;
  bool verbose = false;
  int prev_strand=-1;
  bool first_item = true;


  std::cout << "Parsing JSON file..." << std::endl;
  // Part 1: parse through the JSON file, and make lists called "Scaffolds" and "staples" containing the scaffolds and staples, their strand and index numbers (added custom), and their 3' and 5' neighbors 
  // formatted: Staples: vector<vector<int> >, e.g. Staple[i]={3,27, 3, 26, 3, 28 } would be a staple on strand 3 ind 27 with neighbors at 26 and 28
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



  std::cout << "Finished parsing JSON file.\n";
  //std::cout << "Finding 5' end of scaffold...\n";
  
  // Next step : Find 3' end in list of all scaffold strands and start searching.
  // Cadnano designates the 3' scaffold end as the one whose 3' neighbor is just [-1][-1] 
  int bead_number= 0;
  int scaf_index = Find(ThreePrimeScaffoldEnd[0], ThreePrimeScaffoldEnd[1], Scaffolds);
  int scaffold_nucleotide_number = 0;
  (Scaffolds[scaf_index]).push_back(scaffold_nucleotide_number);
  Scaffolds[scaf_index].push_back(bead_number);

  //std::cout << "Found 5' end, now populating Scaffolds vector...\n";
  while (Scaffolds[scaf_index][4]!=-1 && Scaffolds[scaf_index][5]!=-1){
    //std::cout << "Scaffold number " << scaffold_nucleotide_number << ": " << Scaffolds[scaffold_nucleotide_number][0] << ", " << Scaffolds[scaffold_nucleotide_number][1] << "\n";
    scaffold_nucleotide_number++;
    bead_number = (scaffold_nucleotide_number) / n_nucleotides_per_bead;
    //std::cout << "Bead number " << bead_number << ".\n";
    if (bead_number > (int)((scaffold_nucleotide_number-1)/n_nucleotides_per_bead)){
      connectivity_[bead_number][bead_number-1] += 1;
      connectivity_[bead_number-1][bead_number] += 1;
    }
    scaf_index = Find(Scaffolds[scaf_index][4], Scaffolds[scaf_index][5], Scaffolds);
    Scaffolds[scaf_index].push_back(scaffold_nucleotide_number);
    Scaffolds[scaf_index].push_back(bead_number);

    
    // if (Scaffolds[scaf_index][0] != Scaffolds[scaf_index][4] || Scaffolds[scaf_index][0]!=Scaffolds[scaf_index][2]){
    //  isCrossover[bead_number] = true;
    // }


    //SET STRAND NUMBER
    StrandNumber[bead_number] = Scaffolds[scaf_index][0];
    //////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    /////  STRAND NUMBER HERE
    //////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    
    if (Scaffolds[scaf_index][4]==-1 && Scaffolds[scaf_index][5]==-1 && Scaffolds[scaf_index][2]!=-1){
      if (Scaffolds[scaf_index][4]!=Scaffolds[scaf_index][0]){
	isCrossover[bead_number]=false;
      }
      //std::cout << "found one at " << bead_number <<".\n";
    }
  } //while (Scaffolds!=end) loop done

  size_t n_scaf = bead_number+1;
  for (size_t i = 1; i < n_scaf; ++i){
    if (StrandNumber[i] != StrandNumber[i-1]){
      isCrossover[i]=true;
      isCrossover[i-1]=true;
    }
  }
  if (StrandNumber[0]!=StrandNumber[bead_number] && circular_scaffold){
    isCrossover[0] = true;
    isCrossover[bead_number] = true;
  }

 
  // the next particle in our system will be #n_scaf+1 or (#n_scaf) if 0 is first index.
  int staple_nucleotide_number = scaffold_nucleotide_number+1;

  // Start at the 3' end of each staple that we have identified.
  //std::cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n//identifying connections of " << ThreePrimeStapleEnds.size() << " staples with the scaffold.\n"; //DEBUG 


  // You have to repeat this procedure for each staple in the system separately.
  n_staple_seqs_ = ThreePrimeStapleEnds.size();
  for (size_t i = 0; i < ThreePrimeStapleEnds.size(); i++){
    //std::cout << "Defining staple number " << i << ".\n";
    // This index is the index in our vector of Staple structures containing the three prime staple end.
    int index = Find(ThreePrimeStapleEnds[i][0], ThreePrimeStapleEnds[i][1], Staples);
    // in the structure, record which staple in the simulation will correspond to that part of the design.
    bead_number++;
    stapleNumbers_[bead_number]=i;
    Staples[index].push_back(staple_nucleotide_number);
    int first_nt_number = staple_nucleotide_number;
    Staples[index].push_back(bead_number);
    // If the current staple nucleotide is paired to a scaffold piece, make sure they are designated as connected in the model. 
    int paired_scaf_idx;
    if (index!=-1){
      paired_scaf_idx = Find(Staples[index][0], Staples[index][1], Scaffolds);
    } else {
      paired_scaf_idx = -1;
    }
    if (index!=-1 && paired_scaf_idx!=-1){
      // here Scaffolds[paired_scaf_idx][6] is the in-simulation bead index of the scaffold nucleotide that is connected to this staple.
      SM_[(int)(Scaffolds[paired_scaf_idx][6]/n_nucleotides_per_bead)][bead_number] += 1;
      n_bonds_++;
      //std::cout << "SM[" << (int)(Scaffolds[paired_scaf_idx][6]/n_nucleotides_per_bead) << "][" << bead_number << "] += 1;" << std::endl;
    }
    while(Staples[index][4]!=-1 && Staples[index][5]!=-1){ //As long as current staple nucleotide has a 5' neighbor
      if ((Staples[index][0] != Staples[index][4] && Staples[index][4]!=-1) || (Staples[index][0]!=Staples[index][2] && Staples[index][2]!=-1)){
	//std::cout << "Staple crossover at " << bead_number << "\n";
	isCrossover[bead_number] = true;
      }
      StrandNumber[bead_number] = Staples[index][0];
      // We now know that there is at least one extra staple, so increment staple_nucleotide_number.
      // add the connection between the previous bead and this one.
      //std::cout << "staple_connections[" << bead_number << "][" << bead_number+1 << "] += 1;"<<std::endl;
      staple_nucleotide_number++;

      // If the neighbor accounts for the 11th bead in the series, we need to make a new bead. We also need to add a backbone connection between the previous staple and the new one.
      if ((int)((staple_nucleotide_number-first_nt_number)/n_nucleotides_per_bead) > (int)(((staple_nucleotide_number-1-first_nt_number))/n_nucleotides_per_bead)){
	bead_number++;
	stapleNumbers_[bead_number]=i;
	staple_connections_[bead_number-1][bead_number]+=1;
	staple_connections_[bead_number][bead_number-1]+=1;
	n_bonds_++;
	//std::cout << "staple_connections[" << bead_number-1 << "][" << bead_number << "] = 1;"<<std::endl;
      }
      // get the index of the 5' neighbor staple and set it to the current index.
      index = Find(Staples[index][4], Staples[index][5], Staples);

      // Add the simulation bead # to the design nucleotide
      Staples[index].push_back(staple_nucleotide_number);
      Staples[index].push_back(bead_number);
            
      // Look for scaffold pairing and add bonds as necessary.
      //std::cout << "//Looking for a scaffold at position " << Staples[index][0] <<", " << Staples[index][1] << ".\n";
      paired_scaf_idx = Find(Staples[index][0], Staples[index][1], Scaffolds);
      if (paired_scaf_idx!=-1){
	//std::cout << "debug: paired scaffold index = " << paired_scaf_idx << ", that scaffold bead's bead number = " << Scaffolds[paired_scaf_idx][6] << ".\n"; 
	SM_[(int)(Scaffolds[paired_scaf_idx][7])][bead_number]+=1;
	//std::cout << "SM[" << (int)(Scaffolds[paired_scaf_idx][7]) << "][" << bead_number << "] += 1;"<<std::endl;  
      }	
    }
  }
}
