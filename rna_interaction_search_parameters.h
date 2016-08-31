#ifndef RNA_INTERACTION_SEARCH_PARAMETERS_H
#define RNA_INTERACTION_SEARCH_PARAMETERS_H

#include <string>
#include <cstdio>
#include <iostream>
#include <vector>
using namespace std;

class RnaInteractionSearchParameters{
 private:
  string _db_filename;
  string _input_filename;
  string _output_filename;
  int _hash_size;
  int _repeat_flag;
  int _maximal_span;
  int _min_accessible_length;
  int _max_seed_length;
  double _interaction_energy_threshold;
  double _hybrid_energy_threshold;
  int _drop_out_length_wo_gap;
  int _drop_out_length_w_gap;
  
 public:
  RnaInteractionSearchParameters(){
    _db_filename = "";
    _input_filename = "";
    _hash_size = 0;
    _repeat_flag = 0;
    _maximal_span = 0;
     _min_accessible_length = 0;
    _max_seed_length = 20;
    _hybrid_energy_threshold = -6.5;
    _interaction_energy_threshold = -3.059;
    _drop_out_length_wo_gap = 5;
    _drop_out_length_w_gap = 18;
  }
  void SetParameters(int argc, char* argv[]);
  void SetDbParameters();
  
  string GetDbFilename() const {
    return _db_filename;
  }

  string GetInputFilename() const {
    return _input_filename;
  }

  string GetOutputFilename() const {
    return _output_filename;
  }

  int GetHashSize() const {
    return _hash_size;
  }

  int GetRepeatFlag() const {
    return _repeat_flag;
  }

  int GetMaximalSpan() const {
    return _maximal_span;
  }

  int GetMinAccessibleLength() const {
    return _min_accessible_length;
  }

  int GetMaxSeedLength() const {
    return _max_seed_length;
  }

  double GetInteractionEnergyThreshold() const {
    return _interaction_energy_threshold;
  }

  double GetHybridEnergyThreshold() const {
    return _hybrid_energy_threshold;
  }

  int GetDropOutLengthWoGap() const {
    return _drop_out_length_wo_gap;
  }
  
  int GetDropOutLengthWGap() const {
    return _drop_out_length_w_gap;
  }

  void SetHashSize(int a) {
    _hash_size = a;
  }

  void SetRepeatFlag(int a){
    _repeat_flag = a;
  }

  void SetMaximalSpan(int a) {
    _maximal_span = a;
  }

  void SetMinAccessibleLength(int a) {
    _min_accessible_length = a;
  }
};

#endif
