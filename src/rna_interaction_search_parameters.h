/*
 * rna_interaction_search_parameters.h
 *
 *     Created on: 2016/8/31
 *  Last modified: 2016/11/17
 *         Author: Tsukasa Fukunaga
 */

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
  int _output_style;
  int _hash_size;
  int _repeat_flag;
  int _maximal_span;
  int _min_accessible_length;
  int _max_seed_length;
  double _interaction_energy_threshold;
  double _hybrid_energy_threshold;
  double _final_threshold;
  int _drop_out_length_wo_gap;
  int _drop_out_length_w_gap;
  int _min_helix_length;
 public:
  RnaInteractionSearchParameters(){
    _db_filename = "";
    _input_filename = "";
    _hash_size = 0;
    _repeat_flag = 0;
    _maximal_span = 0;
    _min_accessible_length = 0;
    _output_style = 0;
    _max_seed_length = 20;
    _hybrid_energy_threshold = -6.0;
    _interaction_energy_threshold = -4;
    _final_threshold = -8.0;
    _drop_out_length_wo_gap = 5;
    _drop_out_length_w_gap = 16;
    _min_helix_length = 3;
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

  double GetFinalThreshold() const {
    return _final_threshold;
  }
  
  int GetDropOutLengthWoGap() const {
    return _drop_out_length_wo_gap;
  }
  
  int GetDropOutLengthWGap() const {
    return _drop_out_length_w_gap;
  }

  int GetOutputStyle() const {
    return _output_style;
  }

  int GetMinHelixLength() const {
    return _min_helix_length;
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
