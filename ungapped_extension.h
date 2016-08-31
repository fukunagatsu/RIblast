/*
 * ungapped_extension.h
 *
 *  Created on: 2016/8/31
 *      Author: Tsukasa Fukunaga
 */

#ifndef UNGAPPED_EXTENSION_H
#define UNGAPPED_EXTENSION_H

#include "hit.h"
#include <vector>

using namespace std;

class UngappedExtension {
 public:
  UngappedExtension(int a, int x){
    _min_accessible_length = a;
    _drop_out_score = x;
  }
  void Run(vector<Hit> &candidate, vector<unsigned char> &query_seq, vector<unsigned char> &db_seq, vector<float> &query_accessibility, vector<float> &query_conditional_accessibility, vector<vector<float> > &db_accessibility,vector<vector<float> > &db_conditional_accessibility);

 private:
  int _min_accessible_length;
  int _drop_out_score;
  double LoopEnergy(int type, int type2,int i,int j,int p,int q, vector<unsigned char> &query_seq, vector<unsigned char> &db_seq);
};

#endif
