/*
 * rna_interaction_search.cpp
 *
 *  Created on: 2016/8/31
 *      Author: Tsukasa Fukunaga
 */

#include "rna_interaction_search.h"
#include "fastafile_reader.h"
#include "encoder.h"
#include "sais.h"
#include "raccess.h"
#include "seed_search.h"
#include "ungapped_extension.h"
#include "gapped_extension.h"
#include <fstream>
#include <algorithm>

bool compare( const Hit& left, const Hit& right ) {
  if(left.GetDbSp() != right.GetDbSp()){
    return(left.GetDbSp() < right.GetDbSp());
  }else if(left.GetQSp() != right.GetQSp()){
    return(left.GetQSp() < right.GetQSp());
  }else if(left.GetDbLength() != right.GetDbLength()){
    return(left.GetDbLength() > right.GetDbLength());
  }else{
    return(left.GetQLength() > right.GetQLength());
  }
}

bool base_pair_compare( const BasePair& left, const BasePair& right ) {
  return(left.qpos < right.qpos);
}

void RnaInteractionSearch::Run(const RnaInteractionSearchParameters parameters){
  
  string query_sequence;
  string query_name;
  vector<unsigned char> query_encoded_sequence; query_encoded_sequence.reserve(10);
  vector<int> query_suffix_array;
  vector<float> query_accessibility; query_accessibility.reserve(10);
  vector<float> query_conditional_accessibility; query_conditional_accessibility.reserve(10);
  
  LoadDatabase(parameters.GetDbFilename(), parameters.GetHashSize());
  ReadFastaFile(parameters, query_sequence, query_name);
  CalculateAccessibility(parameters, query_sequence, query_accessibility, query_conditional_accessibility);
  ConstructSuffixArray(parameters, query_sequence, query_encoded_sequence, query_suffix_array);
  vector<Hit> hit_result; hit_result.reserve(10000000);
  
  SearchSeed(parameters,hit_result, query_encoded_sequence, query_suffix_array, query_accessibility, query_conditional_accessibility);

  ExtendWithoutGap(parameters,hit_result, query_encoded_sequence, query_accessibility, query_conditional_accessibility);
  ExtendWithGap(parameters,hit_result, query_encoded_sequence, query_accessibility, query_conditional_accessibility);
  Output(parameters,hit_result,query_name);
}

void RnaInteractionSearch::ReadFastaFile(const RnaInteractionSearchParameters parameters, string &query_sequence, string &query_name){
  FastafileReader fastafile_reader;
  fastafile_reader.ReadFastafile(parameters.GetInputFilename(), query_sequence, query_name);
};

void RnaInteractionSearch::CalculateAccessibility(const RnaInteractionSearchParameters parameters, string &query_sequence, vector<float> &query_accessibility, vector<float> &query_conditional_accessibility){
  FastafileReader fastafile_reader;
  Raccess raccess(parameters.GetMaximalSpan(), parameters.GetMinAccessibleLength());
  raccess.Run(query_sequence, query_accessibility, query_conditional_accessibility);
};

void RnaInteractionSearch::ConstructSuffixArray(const RnaInteractionSearchParameters parameters, string &query_sequence,  vector<unsigned char> &query_encoded_sequence, vector<int> &query_suffix_array){
  Encoder encoder(parameters.GetRepeatFlag());
  encoder.Encode(query_sequence, query_encoded_sequence, 0);
  query_suffix_array.resize(query_encoded_sequence.size());
  sais(&query_encoded_sequence[0], &query_suffix_array[0], query_encoded_sequence.size());
};

void RnaInteractionSearch::SearchSeed(const RnaInteractionSearchParameters parameters, vector<Hit> &hit_result, vector<unsigned char> &query_encoded_sequence, vector<int> &query_suffix_array, vector<float> &query_accessibility, vector<float> &query_conditional_accessibility){
  SeedSearch seed_search(parameters.GetHashSize(), parameters.GetMaxSeedLength(),  parameters.GetMinAccessibleLength(), parameters.GetHybridEnergyThreshold());

  seed_search.Run(query_encoded_sequence, query_suffix_array, _db_seq, _db_suffix_array, _start_hash, _end_hash);
  seed_search.CalcInteractionEnergy(hit_result, query_suffix_array, _db_suffix_array,query_accessibility, query_conditional_accessibility, _db_accessibility, _db_conditional_accessibility, _db_seq_length, _db_seq_start_position);
}

void RnaInteractionSearch::ExtendWithoutGap(const RnaInteractionSearchParameters parameters, vector<Hit> &hit_result, vector<unsigned char> &query_encoded_sequence, vector<float> &query_accessibility, vector<float> &query_conditional_accessibility){
  UngappedExtension ungapped_extension(parameters.GetMinAccessibleLength(), parameters.GetDropOutLengthWoGap());

  ungapped_extension.Run(hit_result, query_encoded_sequence, _db_seq, query_accessibility, query_conditional_accessibility, _db_accessibility, _db_conditional_accessibility);


  double energy_threshold = parameters.GetInteractionEnergyThreshold();
  sort(hit_result.begin(), hit_result.end(), compare);

  for(int i = 0; i<hit_result.size();i++){
    if(hit_result[i].GetEnergy() > energy_threshold){
      hit_result[i].SetFlag();
    }
    if(!hit_result[i].GetFlag()){
      int a_QSp = hit_result[i].GetQSp();
      int a_DbSp = hit_result[i].GetDbSp();
      int a_QEp = a_QSp+hit_result[i].GetQLength()-1;
      int a_DbEp = a_DbSp+hit_result[i].GetDbLength()-1;

      for(int j = i+1; j<hit_result.size();j++){
	if(!hit_result[j].GetFlag()){
	  int b_DbSp = hit_result[j].GetDbSp();
          if(a_DbEp < b_DbSp){
            break;
          }

          int b_QSp = hit_result[j].GetQSp();
          int b_QEp = b_QSp+hit_result[j].GetQLength()-1;
          int b_DbEp = b_DbSp+hit_result[j].GetDbLength()-1;
          if(a_QEp>=b_QEp && a_QSp <= b_QSp && a_DbEp >= b_DbEp){
            hit_result[j].SetFlag();
          }
	}
      }
    }

  }
  hit_result.erase(remove_if(hit_result.begin(), hit_result.end(), CheckFlag()), hit_result.end());
  GetBasePair(hit_result, query_encoded_sequence);
}

void RnaInteractionSearch::ExtendWithGap(const RnaInteractionSearchParameters parameters, vector<Hit> &hit_result, vector<unsigned char> &query_encoded_sequence, vector<float> &query_accessibility, vector<float> &query_conditional_accessibility){
  GappedExtension gapped_extension(parameters.GetMinAccessibleLength(), parameters.GetDropOutLengthWGap());
  gapped_extension.Run(hit_result, query_encoded_sequence, _db_seq, query_accessibility, query_conditional_accessibility, _db_accessibility, _db_conditional_accessibility, _db_seq_length);

  for(int i = 1; i<hit_result.size();i++){
    hit_result[i].SortBasePair();
  }
  sort(hit_result.begin(), hit_result.end(), compare);
  
  for(int i = 1; i<hit_result.size();i++){
    if(SameHitCheckWithGap(hit_result[i], hit_result[i-1])){
        hit_result[i].SetFlag();
    }
  }
  hit_result.erase(remove_if(hit_result.begin(), hit_result.end(), CheckFlag()), hit_result.end());
}

void RnaInteractionSearch::Output(const RnaInteractionSearchParameters parameters, vector<Hit> &hit_result, string q_name){
  ofstream ofs;
  ofs.open(parameters.GetOutputFilename().c_str(),ios::out);
  if (!ofs){
    cout << "Error: can't open output_file:"+parameters.GetOutputFilename()+"." <<endl;
    exit(1);
  }
  for(int i = 0; i< hit_result.size(); i++){
    int id = hit_result[i].GetDbSeqId();
    int length = _db_seq_length[id];
    int seq_start_position = _db_seq_start_position[id];
    ofs << "Query:" << q_name << endl;
    ofs << "Target:" << _db_seq_name[id] << endl;
    ofs << "Energy:" << hit_result[i].GetEnergy() << endl;
    int basepair_length = hit_result[i].GetBasePairLength();
    for(int j = 0; j<basepair_length;j++){
      int dbpos = (length-1) - ( hit_result[i].GetBasePairSecond(j)- seq_start_position);
      ofs << "(" << hit_result[i].GetBasePairFirst(j) << "," << dbpos << ") " ;
    }
    ofs << endl;
    ofs << endl;
  }
  /*if(hit_result.size() != 0){
    double min_energy = INF;
    int k = -1;
    for(int i = 0; i< hit_result.size(); i++){
      if(hit_result[i].GetEnergy() < min_energy){
	min_energy = hit_result[i].GetEnergy();
	k = i;
      }
    }
    int id = hit_result[k].GetDbSeqId();
    int length = _db_seq_length[id];
    int seq_start_position = _db_seq_start_position[id];
    
    int basepair_length = hit_result[k].GetBasePairLength();
    ofs << "Interaction Energy:"<< min_energy << endl;
    for(int i = 0; i<basepair_length;i++){
      int dbpos = (length-1) - ( hit_result[k].GetBasePairSecond(i)- seq_start_position);
      ofs << "(" << hit_result[k].GetBasePairFirst(i) << "," << dbpos << ")" << endl;
    }
  }else{
    ofs << "Interaction Energy:0.0" << endl;
    ofs << "No detection." << endl;
    }*/
  ofs.close();
}


void RnaInteractionSearch::GetBasePair(vector<Hit> &hit_result, vector<unsigned char> &query_encoded_sequence){
  for(int i = 0; i< hit_result.size(); i++){
    int q_start = hit_result[i].GetQSp();
    int db_start = hit_result[i].GetDbSp();
    int length = hit_result[i].GetQLength();
    for(int j = 0; j<length;j++){
      if(BP_pair[query_encoded_sequence[q_start+j]-1][_db_seq[db_start+j]-1] != 0){
	hit_result[i].AddBasePair(q_start+j, db_start+j);
      }
    }
  }
}

bool RnaInteractionSearch::SameHitCheckWithGap(Hit a, Hit b){
  if(a.GetQSp() == b.GetQSp() && a.GetDbSp() == b.GetDbSp() && a.GetQLength() == b.GetQLength() && a.GetDbLength() == b.GetDbLength() && a.GetBasePairLength() == b.GetBasePairLength()){
    int length = a.GetBasePairLength();
    bool flag = true;
    for(int i = 0; i<length; i++){
      if(a.GetBasePairFirst(i) != b.GetBasePairFirst(i) || a.GetBasePairSecond(i) != b.GetBasePairSecond(i)){
	flag = false;
	break;
      }
    }
      
    return(flag);
  }else{
    return(false);
  }
}

void RnaInteractionSearch::LoadDatabase(string db_file_name, int hash_size){
  //load seq file
  ifstream ifs((db_file_name+".seq").c_str(), ios::in | ios::binary);
  if (!ifs){
    cout << "Error: can't open " << db_file_name << ".seq." <<endl;
    exit(1);
  }
  int temp_i = 0;
  ifs.read(reinterpret_cast<char*>(&temp_i), sizeof(int));
  _number_of_db_seq = temp_i;
  vector<int> temp_iv;
  vector<int>::iterator i_it;
  temp_iv.assign(temp_i, 0.0);
  for(i_it=temp_iv.begin();i_it!=temp_iv.end();i_it++){
    ifs.read(reinterpret_cast<char*>(&*i_it),sizeof(int));
  }
  int temp = 0;
  for(int i = 0; i <temp_iv.size(); i++){
    _db_seq_start_position.push_back(temp);
    temp += temp_iv[i]+1;
    _db_seq_length.push_back(temp_iv[i]);
  }
  
  int count = 0;
  vector<unsigned char>::iterator c_it;
  ifs.read(reinterpret_cast<char*>(&count), sizeof(int));
  
  _db_seq.resize(count);
  for(c_it=_db_seq.begin();c_it!=_db_seq.end();c_it++){
    ifs.read(reinterpret_cast<char*>(&*c_it),sizeof(unsigned char));
  }
  ifs.close();
  
  //load acc file
  ifs.open((db_file_name+".acc").c_str(), ios::in | ios::binary);
  if (!ifs){
    cout << "Error: can't open " << db_file_name << ".acc." <<endl;
    exit(1);
  }
  for(int i = 0; i < _number_of_db_seq; i++){
    int a_count = 0;
    vector<float>::iterator it;
    vector <float> temp;
    vector <float> c_temp;
    ifs.read(reinterpret_cast<char*>(&a_count), sizeof(int));
    temp.assign(a_count, 0.0);
    for(it=temp.begin();it!=temp.end();it++){
      ifs.read(reinterpret_cast<char*>(&*it),sizeof(float));
    }
    ifs.read(reinterpret_cast<char*>(&a_count), sizeof(int));
    c_temp.assign(a_count, 0.0);
    for(it=c_temp.begin();it!=c_temp.end();it++){
      ifs.read(reinterpret_cast<char*>(&*it),sizeof(float));
    }
    _db_accessibility.push_back(temp);
    _db_conditional_accessibility.push_back(c_temp);
  }
  ifs.close();

  //load nam file
  ifs.open((db_file_name+".nam").c_str(), ios::in);
  if (!ifs){
    cout << "Error: can't open " << db_file_name << ".nam." <<endl;
    exit(1);
  }
  string str;
  while(getline(ifs, str)){
    _db_seq_name.push_back(str);
  }
  ifs.close();

  //load ind file
  ifs.open((db_file_name+".ind").c_str(), ios::in | ios::binary);
  if (!ifs){
    cout << "Error: can't open " << db_file_name << ".ind." <<endl;
    exit(1);
  }
  count = 0;
  vector<int>::iterator it;
  ifs.read(reinterpret_cast<char*>(&count), sizeof(int));
  
  _db_suffix_array.resize(count);
  for(it=_db_suffix_array.begin();it!=_db_suffix_array.end();it++){
    ifs.read(reinterpret_cast<char*>(&*it),sizeof(int));
  }
  for(int i = 0; i < hash_size; i++){
    for(it=_start_hash[i].begin();it!=_start_hash[i].end();it++){
      ifs.read(reinterpret_cast<char*>(&*it),sizeof(int));
    }
  }
  
  for(int i = 0; i < hash_size; i++){
    for(it=_end_hash[i].begin();it!=_end_hash[i].end();it++){
      ifs.read(reinterpret_cast<char*>(&*it),sizeof(int));
    }
  }
  ifs.close();
}
