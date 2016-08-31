/*
 * gapped_extension.h
 *
 *  Created on: 2016/8/31
 *      Author: Tsukasa Fukunaga
 */

#ifndef GAPPED_EXTENSION_H
#define GAPPED_EXTENSION_H

#include "hit.h"
#include <vector>

using namespace std;

#define MAX_EXTENSION 100000

class Stem{
public:
  Stem(int a, int b, int c){ first = a; second = b; type = c; }
  void set(int a, int b, int c){ first = a; second = b; type = c; }
  int first;  int second; int type;
};

class Cell{
public:
  Cell(int a, int b, double c, double d){
    first = a; second = b; interaction_energy = c; hybrid_energy = d;
  }
  void set(int a, int b, double c, double d){
    first = a; second = b; interaction_energy = c; hybrid_energy = d;
  }
  int first; int second; double hybrid_energy; double interaction_energy;
};

class CheckStemCandidate{
public:
  CheckStemCandidate(int a){ loop_size = a; }
  bool operator()(const Stem &x) const { return(length-x.first-x.second-2 > loop_size);}
  void SetLength(int a){length = a;}
private:
  int loop_size;
  int length;
};

class GappedExtension {
 public:
  GappedExtension(int a, int x){
    _min_accessible_length = a;
    _drop_out_score = x;
  }
  void Run(vector<Hit> &candidate, vector<unsigned char> &query_seq, vector<unsigned char> &db_seq, vector<float> &query_accessibility, vector<float> &query_conditional_accessibility, vector<vector<float> > &db_accessibility,vector<vector<float> > &db_conditional_accessibility, vector<int> &db_seq_length);

 private:
  int _min_accessible_length;
  int _drop_out_score;

  int GetChar(vector<unsigned char> &seq, int i);
  void traceback(Hit &candidate, vector<vector<Cell> > &matrix_c, int i,int i_start, int j, int j_start, bool flag);
  void extension(Hit &candidate, vector<unsigned char> &query_seq, vector<unsigned char> &db_seq, vector<float> &query_accessibility, vector<float> &query_conditional_accessibility, vector<vector<float> > &db_accessibility,vector<vector<float> > &db_conditional_accessibility, bool flag);
  double LoopEnergy(int type, int type2,int i,int j,int p,int q, vector<unsigned char> &query_seq, vector<unsigned char> &db_seq);
  double CalcDangleEnergy(int q_pos, int db_pos, int flag, vector<unsigned char> &query_seq, vector<unsigned char> &db_seq, int dbseq_length);
};

#endif
