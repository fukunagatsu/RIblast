/*
 * rna_interaction_search.cpp
 *
 *     Created on: 2016/8/31
 *  Last modified: 2016/11/21
 *         Author: Tsukasa Fukunaga
 */

#ifndef RNA_INTERACTION_SEARCH_H
#define RNA_INTERACTION_SEARCH_H

#include "rna_interaction_search_parameters.h"
#include "hit.h"
#include <vector>
#include <math.h>

class CheckFlag{
public:
  bool operator()(const Hit &a) const { return(a.GetFlag()); }
};


class RnaInteractionSearch {
 public:
  RnaInteractionSearch(int hash_size) {
    int nucleotide_size = 4;
    _start_hash.reserve(hash_size);
    _end_hash.reserve(hash_size);
    for(int i = 0; i < hash_size; i++){
      vector <int> temp_vector; temp_vector.resize((int)pow(nucleotide_size,i+1), 0);
      _start_hash.push_back(temp_vector);
      vector <int> temp2_vector; temp_vector.resize((int)pow(nucleotide_size,i+1), 0);
      _end_hash.push_back(temp_vector);
    }
    _coefficient_a = 0.0;
    _coefficient_b = 0.0;
    _eta = 0.0;
  }
  void Run(const RnaInteractionSearchParameters parameters);
 private:
  void ReadFastaFile(const RnaInteractionSearchParameters parameters, vector<string> &sequences, vector<string> &names);
  void CalculateAccessibility(const RnaInteractionSearchParameters parameters, string &query_sequence, vector<float> &query_accessibility, vector<float> &query_conditional_accessibility);
  void ConstructSuffixArray(const RnaInteractionSearchParameters parameters, string &query_sequence,  vector<unsigned char> &query_encoded_sequence, vector<int> &query_suffix_array);
  void SearchSeed(const RnaInteractionSearchParameters parameters, vector<Hit> &hit_result, vector<unsigned char> &query_encoded_sequence, vector<int> &query_suffix_array, vector<float> &query_accessibility, vector<float> &query_conditional_accessibility);
  void ExtendWithoutGap(const RnaInteractionSearchParameters parameters, vector<Hit> &hit_result, vector<unsigned char> &query_encoded_sequence, vector<float> &query_accessibility, vector<float> &query_conditional_accessibility);
  void ExtendWithGap(const RnaInteractionSearchParameters parameters, vector<Hit> &hit_result, vector<unsigned char> &query_encoded_sequence, vector<float> &query_accessibility, vector<float> &query_conditional_accessibility);
  double CalcPvalue(int q_length, int db_length, double energy);
  void CalcPvalueParameter(double w);
  void GetBasePair(vector<Hit> &hit_result, vector<unsigned char> &query_encoded_sequence);
  void Output(const RnaInteractionSearchParameters parameters, vector<Hit> &hit_result, string q_name, int flag, int *count, int q_length);
  void LoadDatabase(string db_file_name, int hash_size);
  void CheckRedundancy(vector<Hit> &hit_result, double energy_threshold);
  bool SameHitCheckWithGap(Hit a, Hit b);
  vector<vector<float> > _db_accessibility;
  vector<vector<float> > _db_conditional_accessibility;
  vector<string> _db_seq_name;
  vector<int> _db_suffix_array;
  vector<unsigned char> _db_seq;
  vector<vector<int> > _start_hash;
  vector<vector<int> > _end_hash;
  vector<int> _db_seq_start_position;
  vector<int> _db_seq_length;
  int _number_of_db_seq;
  double _coefficient_a;
  double _coefficient_b;
  double _eta;
};

#endif
