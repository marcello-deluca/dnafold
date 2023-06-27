#ifndef IMPORTER_FUNCTION_H
#define IMPORTER_FUNCTION_H
#include <string>
void load_file(char * inputfilename_, std::vector<std::vector<int> > & SM_, std::vector<std::vector<int> > & staple_connections_, std::vector<int> & stapleNumbers_, std::vector<std::vector<int> > & connectivity_, size_t & n_bonds_, size_t & n_staple_seqs_, std::vector<bool> & isCrossover, std::vector<int> & StrandNumber, double n_nucleotides_per_beado); 
#endif
