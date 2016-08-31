#ifndef DB_CONSTRUCTION_H
#define DB_CONSTRUCTION_H

#include "db_construction_parameters.h"
#include <vector>

class DbConstruction {
 public:
  DbConstruction() {}
  void Run(const DbConstructionParameters parameters);
 private:
  void ReadFastaFile(const DbConstructionParameters parameters,  vector<string> &sequences);
  void CalculateAccessibility(const DbConstructionParameters parameters, vector<string> &sequences);
  void ConstructSuffixArray(const DbConstructionParameters parameters, vector<string> &sequences,  vector<unsigned char> &encoded_sequences, vector<int> &suffix_array);
  void ConstructHashForShortSubstring(const DbConstructionParameters parameters,  vector<unsigned char> &encoded_sequences, vector<int> &suffix_array, vector<vector <int> > &start_hash, vector<vector <int> > &end_hash);
  void Search(vector<unsigned char> &encoded_sequences, vector<int> &suffix_array, int* start, int* end, unsigned char c, int offset);
  void SaveBasicInformation(DbConstructionParameters parameters, vector<string> &names);
  void SaveIndexData(string file_name, vector<int> &suffix_array, vector<vector <int> > &start_hash, vector<vector <int> > &end_hash);
  void SaveSequenceData(string file_name,  vector<string> &sequences, vector<unsigned char> &encoded_sequences);
};

#endif
