#ifndef HIT_H
#define HIT_H

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;

class BasePair{
public:
  BasePair(int a, int b){ qpos = a; dbpos = b;}
  int qpos;  int dbpos;
};

class Hit{
 public:
  Hit(int q_sp, int db_sp, int length, double energy){
    _dbseq_id = -1;
    _dbseq_id_start = -1;
    _q_sp = q_sp;
    _db_sp = db_sp;
    _q_length = length;
    _db_length = length;
    _energy = energy;
    _flag = false;
    _base_pair.reserve(10);
  }

  int GetDbSeqId() const{
    return _dbseq_id;
  }

  static bool base_pair_compare( const BasePair& left, const BasePair& right ) {
    return(left.qpos < right.qpos);
  }

  int GetDbSeqIdStart() const{
    return _dbseq_id_start;
  }

  int GetQSp() const{
    return _q_sp;
  }

  int GetDbSp() const{
    return _db_sp;
  }

  unsigned short int GetQLength() const{
    return _q_length;
  }

  unsigned short int GetDbLength() const{
    return _db_length;
  }

  int GetBasePairLength() const{
    return _base_pair.size();
  }

  double GetEnergy() const{
    return _energy;
  }
  
  bool GetFlag() const{
    return(_flag);
  }
  
  int GetBasePairFirst(int i) const{
    return(_base_pair[i].qpos);
  }
  
  int GetBasePairSecond(int i) const{
    return(_base_pair[i].dbpos);
  }

  void SetEnergy(double a){
    _energy = a;
  }

  void SetDbSeqId(int a){
    _dbseq_id = a;
  }

  void SetDbSeqIdStart(int a){
    _dbseq_id_start = a;
  }

  void SetFlag(){
    _flag = true;
  }

  void SetQSp(int a){
    _q_sp = a;
  }

  void SetDbSp(int a){
    _db_sp = a;
  }

  void SetQLength(int a){
    _q_length = a;
  }

  void SetDbLength(int a){
    _db_length = a;
  }

  void AddBasePair(int a, int b){
    _base_pair.push_back(BasePair(a,b));
  }
  
  void SortBasePair(){
    sort(_base_pair.begin(), _base_pair.end(), base_pair_compare);
  }
  
 private:
  int _dbseq_id;
  int _dbseq_id_start;
  int _q_sp;
  int _db_sp;
  int _q_length;
  int _db_length;
  double _energy;
  bool _flag;
  vector<BasePair> _base_pair;
};

class Hit_candidate{
 public:
  Hit_candidate(int sp_q_sa, int ep_q_sa, int sp_db_sa, int ep_db_sa, int length, double energy){
    _sp_q_sa = sp_q_sa;
    _ep_q_sa = ep_q_sa;
    _sp_db_sa = sp_db_sa;
    _ep_db_sa = ep_db_sa;
    _length = length;
    _energy = energy;
  }

  double GetEnergy() const{
    return _energy;
  }

  int GetSpQSa() const{
    return _sp_q_sa;
  }

  int GetEpQSa() const{
    return _ep_q_sa;
  }

  int GetSpDbSa() const{
    return _sp_db_sa;
  }

  int GetEpDbSa() const{
    return _ep_db_sa;
  }

  int GetLength() const{
    return _length;
  }
  
 private:
  int _sp_q_sa;
  int _ep_q_sa;
  int _sp_db_sa;
  int _ep_db_sa;
  int _length;
  double _energy;
};

#endif
