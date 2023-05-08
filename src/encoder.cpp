/*
 * encoder.cpp
 *
 *  Created on: 2016/8/31
 *      Author: Tsukasa Fukunaga
 */

#include "encoder.h"

void Encoder::Encode(vector<string> &sequences, vector<unsigned char> &encoded_sequences, int r){
  for(int i = 0; i < sequences.size(); i++){
    if(r == 0){
      for(int j = 0; j<sequences[i].size();j++){
	encoded_sequences.push_back(_code_table[static_cast<int>(sequences[i][j])]);  
      }
    }else{
      for(int j = sequences[i].size()-1; j>= 0; j--){
	encoded_sequences.push_back(_code_table[static_cast<int>(sequences[i][j])]);  
      }
    }
    encoded_sequences.push_back(_sentinel_character);
  }
}

void Encoder::Encode(string &sequence, vector<unsigned char> &encoded_sequence, int r){
  if(r == 0){
    for(int i = 0; i<sequence.size();i++){
      encoded_sequence.push_back(_code_table[static_cast<int>(sequence[i])]);  
    }
  }else{
    for(int i = sequence.size()-1; i>= 0; i--){
      encoded_sequence.push_back(_code_table[static_cast<int>(sequence[i])]);  
    }
  }
  encoded_sequence.push_back(_sentinel_character);
}
