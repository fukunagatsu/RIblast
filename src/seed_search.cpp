/*
 * seed_search.cpp
 *
 *     Created on: 2016/8/31
 *  Last modified: 2016/11/21
 *         Author: Tsukasa Fukunaga
 */

#include "seed_search.h"
#include <math.h>
#include <iostream>
          
void SeedSearch::Run(vector<unsigned char> &query_seq, vector<int> &query_suffix_array, vector<unsigned char> &db_seq, vector<int> &db_suffix_array,  vector<vector<int> > &start_hash,  vector<vector<int> > &end_hash){
  vector<int> db_seed; db_seed.reserve(_max_seed_length);
  vector<int> q_seed; q_seed.reserve(_max_seed_length);
  
  SeedSearchCore(query_seq, query_suffix_array, db_seq, db_suffix_array, start_hash, end_hash, db_seed, q_seed, 0, query_suffix_array.size()-1, 0, db_suffix_array.size()-1, 0.0, 0);
  return;
}

void SeedSearch::CalcInteractionEnergy(vector<Hit> &hit_result, vector<int> &query_suffix_array, vector<int> &db_suffix_array, vector<float> &query_accessibility, vector<float> &query_conditional_accessibility, vector<vector<float> > &db_accessibility,vector<vector<float> > &db_conditional_accessibility, vector<int> &db_seq_length, vector<int> &db_seq_start_position){
  for(int i = 0;i < _hit_candidate_result.size(); i++){
    int sp_q_sa = _hit_candidate_result[i].GetSpQSa();
    int ep_q_sa = _hit_candidate_result[i].GetEpQSa();
    int sp_db_sa = _hit_candidate_result[i].GetSpDbSa();
    int ep_db_sa = _hit_candidate_result[i].GetEpDbSa();
    int length = _hit_candidate_result[i].GetLength();
    double temp_score = _hit_candidate_result[i].GetEnergy();

    vector<int> q_sp_vector; q_sp_vector.reserve(10000);
    vector<double> qa_vector; qa_vector.reserve(10000);
    for(int j = sp_q_sa; j <=ep_q_sa ; j++){
      q_sp_vector.push_back(query_suffix_array[j]);
      qa_vector.push_back(CalcAccessibility(query_accessibility, query_conditional_accessibility, query_suffix_array[j], length));
    }

    for(int k = sp_db_sa; k <=ep_db_sa ; k++){
      int db_sp = db_suffix_array[k];
      int dbseq_start = 0; int dbseq_id = 0;
      GetSeqIdAndStart(db_seq_length, db_seq_start_position, &dbseq_id, &dbseq_start, db_sp, length);
      double dba = CalcAccessibility(db_accessibility[dbseq_id], db_conditional_accessibility[dbseq_id], dbseq_start, length);
      
      for(int j = sp_q_sa; j <= ep_q_sa; j++){
	double interaction_energy = qa_vector[j-sp_q_sa] + dba + temp_score;
	
	if(interaction_energy < 0){
	  Hit temp_hit( q_sp_vector[j-sp_q_sa], db_sp, length, qa_vector[j-sp_q_sa] + dba, temp_score);
	  temp_hit.SetDbSeqId(dbseq_id);
	  temp_hit.SetDbSeqIdStart(dbseq_start);
	  try{
	    hit_result.push_back(temp_hit);
	  }catch(bad_alloc){
	    cerr << "The memory capacity is too small to run RIblast for this dataset" << endl;  
	    exit(1);
	  }
	}
      }
    } 
  }
  return;
}

void SeedSearch::GetSeqIdAndStart(vector<int> &db_seq_length, vector<int> &db_seq_start_position, int* seq_id, int* start, int sp, int length){
  int s = 0;
  int e = db_seq_start_position.size()-1;
  int m = (s+e)/2;

  while(true){
    if(e-s == 0){
      *seq_id  = s; *start = db_seq_length[s] - (sp - db_seq_start_position[s]) - length;
      break;
    }else if(e-s==1){
      if(sp >= db_seq_start_position[s] && sp < db_seq_start_position[s+1]){
	*seq_id  = s; *start = db_seq_length[s] - (sp - db_seq_start_position[s]) - length;
	break;
      }else{
	*seq_id  = e; *start = db_seq_length[e] - (sp - db_seq_start_position[e]) - length;
	break;
      }
    }else{    
      if(sp >= db_seq_start_position[m] && sp < db_seq_start_position[m+1]){
	*seq_id  = m;
	*start = db_seq_length[m] - (sp - db_seq_start_position[m]) - length;
	break;
      }else{
	if(sp >= db_seq_start_position[m]){
	  s = m;
	  m = (s+e)/2;
	}else{
	  e = m;
	  m = (s+e)/2;
	}
      }
    }
  }
  
  return;
}

double SeedSearch::CalcAccessibility(vector<float> &accessibility, vector<float> &conditional_accessibility, int sp, int length){
  double temp = accessibility[sp];
  for(int i =_min_accessible_length; i < length; i++){
    temp += conditional_accessibility[sp+i];
  }
  return(temp);
}

void SeedSearch::SeedSearchCore(vector<unsigned char> &query_seq, vector<int> &query_suffix_array, vector<unsigned char> &db_seq, vector<int> &db_suffix_array,  vector<vector<int> > &start_hash,  vector<vector<int> > &end_hash, vector<int> &db_seed, vector<int> &q_seed, int sp_q, int ep_q, int sp_db, int ep_db, double score, int length){

  if(length < _max_seed_length){
    vector<int> sp_q_result; sp_q_result.reserve(_pair_size);
    vector<int> ep_q_result; ep_q_result.reserve(_pair_size);
    vector<int> sp_db_result; sp_db_result.reserve(_pair_size);
    vector<int> ep_db_result; ep_db_result.reserve(_pair_size);

    for(int i = 0; i< _pair_size; i++){
      int start = sp_q;
      int end = ep_q;
      int nucleotide_size = 4;
      SeedSearchNextCharacter(query_seq, query_suffix_array, &start, &end, _stem_pair[i][0], length);
      sp_q_result.push_back(start);
      ep_q_result.push_back(end);
      start = sp_db;
      end = ep_db;
      if(length+1 > _hash_size){
	SeedSearchNextCharacter(db_seq, db_suffix_array, &start, &end, _stem_pair[i][1], length);
      }else{
	int temp = _stem_pair[i][1]-2;
	for(int j = 0; j <length;j++){
	  temp += (int)pow(nucleotide_size,length-j) * (db_seed[j] - 2);
	}
	start = start_hash[length][temp];
	end = end_hash[length][temp];
      }
      sp_db_result.push_back(start);
      ep_db_result.push_back(end);
    }
    
    for(int i = 0; i< _pair_size; i++){
      if(sp_q_result[i] <= ep_q_result[i] && sp_db_result[i] <= ep_db_result[i]){
	double temp_score = 0.0;
	
	if(length > 0){
	  int type = BP_pair[q_seed[length-1]-1][db_seed[length-1]-1];
	  int type2 = BP_pair[_stem_pair[i][0]-1][_stem_pair[i][1]-1];
	  type2 = rtype[type2];
	  temp_score = score + ((double)stack37[type][type2])/100;
	}
	
	if(temp_score < _hybrid_energy_threshold && length+1 >= _min_accessible_length){
	  Hit_candidate temp_hit_candidate(sp_q_result[i], ep_q_result[i], sp_db_result[i], ep_db_result[i], length+1, temp_score);
	  _hit_candidate_result.push_back(temp_hit_candidate);
	}else{
	  q_seed.push_back(_stem_pair[i][0]);
	  db_seed.push_back(_stem_pair[i][1]);
	  SeedSearchCore(query_seq, query_suffix_array, db_seq, db_suffix_array, start_hash, end_hash, db_seed, q_seed, sp_q_result[i], ep_q_result[i], sp_db_result[i], ep_db_result[i], temp_score, length+1);
	}
      }
    }
  }
  q_seed.pop_back();
  db_seed.pop_back();
  return;
}

void SeedSearch::SeedSearchNextCharacter(vector<unsigned char> &encoded_sequences, vector<int> &suffix_array, int* start, int* end, unsigned char c, int offset){
  int s = *start;
  int e = *end;
  int m;

  if (suffix_array[s] + offset >= encoded_sequences.size()) {
    ++(*start);
  }
  
  if(s>e){
    *start = 1;
    *end = 0;
    return;
  }else if (s == e) {
    if (encoded_sequences[suffix_array[s]+offset] == c) {
      return;
    } else {
      *start = 1;
      *end = 0;
      return;
    }
  }
  
  if (encoded_sequences[suffix_array[s]+offset] != c) {
    while (s < e - 1) {
      m = (s + e) / 2;
      if (encoded_sequences[suffix_array[m]+offset] < c) {
        s = m;
      } else {
        e = m;
      }
    }
    if (encoded_sequences[suffix_array[e]+offset] != c) {
      *start = 1;
      *end = 0;
      return;
    }
    *start = e;
    s = e;
    e = *end;
  }

  if (encoded_sequences[suffix_array[e]+offset] != c) {
    while (s < e - 1) {
      m = (s + e) / 2;
      if (encoded_sequences[suffix_array[m]+offset] > c) {
        e = m;
      } else {
        s = m;
      }
    }
    if (encoded_sequences[suffix_array[s]+offset] != c) {
      *start = 1;
      *end = 0;
      return;
    }
    *end = s;
  }
  return;
}
