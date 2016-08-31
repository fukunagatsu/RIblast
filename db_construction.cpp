/*
 * db_construction.cpp
 *
 *  Created on: 2016/8/31
 *      Author: Tsukasa Fukunaga
 */

#include "db_construction.h"
#include "fastafile_reader.h"
#include "encoder.h"
#include "sais.h"
#include "raccess.h"
#include <math.h>
#include <time.h> 
          
void DbConstruction::Run(const DbConstructionParameters parameters){
  vector<string> sequences; sequences.reserve(10);
  vector<unsigned char> encoded_sequences;
  vector<int> suffix_array;
  vector<vector <int> > start_hash; start_hash.reserve(parameters.GetHashSize());
  vector<vector <int> > end_hash; start_hash.reserve(parameters.GetHashSize());

  ReadFastaFile(parameters, sequences);
  CalculateAccessibility(parameters,sequences);
  ConstructSuffixArray(parameters, sequences,encoded_sequences,suffix_array);
  ConstructHashForShortSubstring(parameters, encoded_sequences, suffix_array, start_hash, end_hash);

  SaveIndexData(parameters.GetDbFilename(), suffix_array, start_hash, end_hash);
  SaveSequenceData(parameters.GetDbFilename(), sequences, encoded_sequences);
  
}

void DbConstruction::ReadFastaFile(const DbConstructionParameters parameters,  vector<string> &sequences){
  vector<string> names; names.reserve(10);
  FastafileReader fastafile_reader;
  fastafile_reader.ReadFastafile(parameters.GetInputFilename(), sequences, names);
  SaveBasicInformation(parameters, names);
}

void DbConstruction::CalculateAccessibility(const DbConstructionParameters parameters, vector<string> &sequences){
  Raccess raccess(parameters.GetDbFilename(), parameters.GetMaximalSpan(), parameters.GetMinAccessibleLength());
  for (int i = 0; i < sequences.size() ; i++){
    raccess.Run(sequences[i]);
  }
}

void DbConstruction::ConstructSuffixArray(const DbConstructionParameters parameters, vector<string> &sequences,  vector<unsigned char> &encoded_sequences, vector<int> &suffix_array){
  Encoder encoder(parameters.GetRepeatFlag());
  encoder.Encode(sequences, encoded_sequences, 1);
  suffix_array.resize(encoded_sequences.size());
  sais(&encoded_sequences[0], &suffix_array[0], encoded_sequences.size());
};

void DbConstruction::ConstructHashForShortSubstring(const DbConstructionParameters parameters,  vector<unsigned char> &encoded_sequences, vector<int> &suffix_array, vector<vector <int> > &start_hash, vector<vector <int> > &end_hash){
  int nucleotide_size = 4;
  
  for(int i = 0; i < parameters.GetHashSize(); i++){
    vector<int> temp_vector; temp_vector.reserve((int)pow(nucleotide_size,i+1));
    start_hash.push_back(temp_vector);
    vector<int> temp_vector2; temp_vector2.reserve((int)pow(nucleotide_size,i+1));
    end_hash.push_back(temp_vector2);
  }

  int start = 0;  int end = 0;
  for(int i = 0; i < parameters.GetHashSize(); i++){
    for(int j  =0; j < (int)pow(nucleotide_size,i+1); j++){
      unsigned char c = (j%4)+2;
      if(i == 0){
	start = 0;
	end = suffix_array.size()-1;
      }else{
	start = start_hash[i-1][j/4];
	end = end_hash[i-1][j/4];
      }
      Search(encoded_sequences, suffix_array, &start, &end, c, i);
      start_hash[i].push_back(start);
      end_hash[i].push_back(end);
    }
  }
}

void DbConstruction::SaveSequenceData(string file_name, vector<string> &sequences, vector<unsigned char> &encoded_sequences){
  ofstream of((file_name+".seq").c_str(), ios::out | ios::binary);
  int number_of_seq = sequences.size();
  of.write(reinterpret_cast<const char*>(&number_of_seq), sizeof(int));
  for(int i = 0; i < number_of_seq; i++){
    int seq_size = sequences[i].size();
    of.write(reinterpret_cast<const char*>(&seq_size), sizeof(int));
  }
  int count = encoded_sequences.size();
  vector<unsigned char>::iterator it;
  of.write(reinterpret_cast<const char*>(&count), sizeof(int));
  for(it = encoded_sequences.begin(); it != encoded_sequences.end(); it++){
    of.write(reinterpret_cast<const char*>(&*it), sizeof(unsigned char));
  }
  of.close();
}

void  DbConstruction::SaveIndexData(string file_name, vector<int> &suffix_array, vector<vector <int> > &start_hash, vector<vector <int> > &end_hash){
  ofstream of((file_name+".ind").c_str(), ios::out | ios::binary);
  int count = suffix_array.size();
  vector<int>::iterator it;
  of.write(reinterpret_cast<const char*>(&count), sizeof(int));
  for(it = suffix_array.begin(); it != suffix_array.end(); it++){
    of.write(reinterpret_cast<const char*>(&*it), sizeof(int));
  }
  for(int i = 0; i< start_hash.size();i++){
    for(it = start_hash[i].begin(); it != start_hash[i].end(); it++){
      of.write(reinterpret_cast<const char*>(&*it), sizeof(int));
    }
  }
  for(int i = 0; i< end_hash.size();i++){
    for(it = end_hash[i].begin(); it != end_hash[i].end(); it++){
      of.write(reinterpret_cast<const char*>(&*it), sizeof(int));
    }
  }
  of.close();
}

void DbConstruction::SaveBasicInformation(DbConstructionParameters parameters, vector<string> &names){
  ofstream of((parameters.GetDbFilename()+".bas").c_str(), ios::out | ios::binary);
  int hash_size = parameters.GetHashSize();
  int repeat_flag = parameters.GetRepeatFlag();
  int maximal_span = parameters.GetMaximalSpan();
  int min_accessible_length = parameters.GetMinAccessibleLength();
  of.write(reinterpret_cast<const char*>(&hash_size), sizeof(int));
  of.write(reinterpret_cast<const char*>(&repeat_flag), sizeof(int));
  of.write(reinterpret_cast<const char*>(&maximal_span), sizeof(int));
  of.write(reinterpret_cast<const char*>(&min_accessible_length), sizeof(int));

  of.close();
  of.open((parameters.GetDbFilename()+".nam").c_str(), ios::out);
  for(int i = 0; i < names.size(); i++){
    of << names[i] << endl;
  }

  of.close();
}

void DbConstruction::Search(vector<unsigned char> &encoded_sequences, vector<int> &suffix_array, int* start, int* end, unsigned char c, int offset){
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
