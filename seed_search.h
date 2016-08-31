#ifndef SEED_SEARCH_H
#define SEED_SEARCH_H

#include <vector>
#include "energy_par.h"
#include "hit.h"

using namespace std;

class SeedSearch {
 public:
  SeedSearch(int h, int l, int mal, double he){
    _max_seed_length = l;
    _pair_size = 6;
    _hash_size = h;
    _min_accessible_length = mal;
    _hybrid_energy_threshold = he;
    
    _stem_pair.resize(_pair_size);
    _stem_pair[0].push_back(3); _stem_pair[0].push_back(4);
    _stem_pair[1].push_back(4); _stem_pair[1].push_back(3);
    _stem_pair[2].push_back(4); _stem_pair[2].push_back(5);
    _stem_pair[3].push_back(5); _stem_pair[3].push_back(4);
    _stem_pair[4].push_back(2); _stem_pair[4].push_back(5);
    _stem_pair[5].push_back(5); _stem_pair[5].push_back(2);
    _hit_candidate_result.reserve(5000);
  }
  void Run(vector<unsigned char> &query_seq, vector<int> &query_suffix_array, vector<unsigned char> &db_seq, vector<int> &db_suffix_array,  vector<vector<int> > &_start_hash,  vector<vector<int> > &_end_hash);
  void CalcInteractionEnergy(vector<Hit> &hit_result, vector<int> &query_suffix_array, vector<int> &db_suffix_array, vector<float> &query_accessibility, vector<float> &query_conditional_accessibility, vector<vector<float> > &db_accessibility,vector<vector<float> > &db_conditional_accessibility, vector<int> &_db_seq_length, vector<int> &_db_seq_start_position);
  
 private:
  vector<vector<int> > _stem_pair;
  int _pair_size;
  int _max_seed_length;
  int _hash_size;
  int _min_accessible_length;
  double _hybrid_energy_threshold;
  vector<Hit_candidate> _hit_candidate_result;

  void GetSeqIdAndStart(vector<int> &db_seq_length, vector<int> &db_seq_start_position, int* seq_id, int* start, int sp, int length);
  void SeedSearchCore(vector<unsigned char> &query_seq, vector<int> &query_suffix_array, vector<unsigned char> &db_seq, vector<int> &db_suffix_array,  vector<vector<int> > &_start_hash,  vector<vector<int> > &_end_hash,vector<int> &db_seed, vector<int> &q_seed, int sp_q, int ep_q, int sp_db, int ep_db, double score, int length);
  void SeedSearchNextCharacter(vector<unsigned char> &encoded_sequences, vector<int> &suffix_array, int* start, int* end, unsigned char c, int offset);
  double CalcAccessibility(vector<float> accessibility, vector<float> conditional_accessibility, int sp, int length);
};

#endif
