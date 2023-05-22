/*
 * encoder.h
 *
 *  Created on: 2016/8/31
 *      Author: Tsukasa Fukunaga
 */

#ifndef ENCODER_H
#define ENCODER_H

#include <limits.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <iostream>

using namespace std;

class Encoder {
 public:
  Encoder(int flag) {
    _sentinel_character = 0;
    _unknown_character = 1; 
    if(flag == 0){
      for (int i = 0; i < UCHAR_MAX; ++i) {
	_code_table[i] = _unknown_character;
      }
      
      _code_table[static_cast<int>('A')] = 2;
      _code_table[static_cast<int>('C')] = 3;
      _code_table[static_cast<int>('G')] = 4;
      _code_table[static_cast<int>('T')] = 5;
      _code_table[static_cast<int>('U')] = 5;
      
    }else if(flag == 1){
      for (int i = 0; i < UCHAR_MAX; ++i) {
	_code_table[i] = _unknown_character;
      }
      
      _code_table[static_cast<int>('A')] = 2;
      _code_table[static_cast<int>('C')] = 3;
      _code_table[static_cast<int>('G')] = 4;
      _code_table[static_cast<int>('T')] = 5;
      _code_table[static_cast<int>('U')] = 5;
      _code_table[static_cast<int>('a')] = 6;
      _code_table[static_cast<int>('c')] = 7;
      _code_table[static_cast<int>('g')] = 8;
      _code_table[static_cast<int>('t')] = 9;
      _code_table[static_cast<int>('u')] = 9;
      
    }else if(flag == 2){
      for (int i = 0; i < UCHAR_MAX; ++i) {
	_code_table[i] = _unknown_character;
      }
      _code_table[static_cast<int>('A')] = 2;
      _code_table[static_cast<int>('C')] = 3;
      _code_table[static_cast<int>('G')] = 4;
      _code_table[static_cast<int>('T')] = 5;
      _code_table[static_cast<int>('U')] = 5;
      _code_table[static_cast<int>('a')] = 2;
      _code_table[static_cast<int>('c')] = 3;
      _code_table[static_cast<int>('g')] = 4;
      _code_table[static_cast<int>('t')] = 5;
      _code_table[static_cast<int>('u')] = 5;
      
    }else{
      cerr << "Error: -r option must be 0, 1, or 2." << endl;
      exit(1);
    }
  }
  void Encode(vector<string> &sequences, vector<unsigned char> &encoded_sequences, int r);
  void Encode(string &sequence, vector<unsigned char> &encoded_sequences, int r);
 private:
  unsigned char _code_table[UCHAR_MAX];
  unsigned char _sentinel_character;
  unsigned char _unknown_character;
};

#endif
