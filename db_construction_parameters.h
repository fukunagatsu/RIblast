/*
 * db_construction_parameters.h
 *
 *  Created on: 2016/8/31
 *      Author: Tsukasa Fukunaga
 */

#ifndef DB_CONSTRUCTION_PARAMETERS_H
#define DB_CONSTRUCTION_PARAMETERS_H

#include <string>
#include <cstdio>
#include <iostream>
using namespace std;

class DbConstructionParameters{
 private:
  string _db_filename;
  string _input_filename;
  int _hash_size;
  int _repeat_flag;
  int _maximal_span;
  int _min_accessible_length;
  double _accessibility_threshold;
  
 public:
  DbConstructionParameters(){
    _db_filename = "";
    _input_filename = "";
    _hash_size = 8;
    _repeat_flag = 0;
    _maximal_span = 70;
    _min_accessible_length = 5;
    _accessibility_threshold = 4.5;
  }
  void SetParameters(int argc, char* argv[]);
  
  string GetDbFilename() const {
    return _db_filename;
  }

  string GetInputFilename() const {
    return _input_filename;
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

  double GetAccessibilityThreshold() const{
    return _accessibility_threshold;
  }
};

#endif
