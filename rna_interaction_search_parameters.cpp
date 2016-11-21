/*
 * rna_interaction_search_parameters.cpp
 *
 *     Created on: 2016/8/31
 *  Last modified: 2016/11/17
 *         Author: Tsukasa Fukunaga
 */

#include "rna_interaction_search_parameters.h"
#include <getopt.h>
#include <stdlib.h>
#include <fstream>

void RnaInteractionSearchParameters::SetParameters(int argc, char* argv[]) {
  int c;
  extern char *optarg;
  while ((c = getopt(argc, argv, "i:o:d:l:e:y:x:f:g:s:")) != -1) {
    switch (c) {
    case 'i':
      _input_filename = optarg;
      break;

    case 'o':
      _output_filename = optarg;
      break;
      
    case 'd':
      _db_filename = optarg;
      break;

    case 'l':
      _max_seed_length = atoi(optarg);
      break;

    case 'e':
      _hybrid_energy_threshold = atof(optarg);
      break;
   
    case 'f':
      _interaction_energy_threshold = atof(optarg);
      break;

    case 'g':
      _final_threshold = atof(optarg);
      break;

    case 's':
      _output_style = atoi(optarg);
      break;

    case 'x':
      _drop_out_length_w_gap = atoi(optarg);
      break;

    case 'y':
      _drop_out_length_wo_gap = atoi(optarg);
      break;

    default:
      cerr << "Error: The argument is invalid command." << endl;
      exit(1);
    }
  }
}

void RnaInteractionSearchParameters::SetDbParameters(){
  ifstream ifs((GetDbFilename()+".bas").c_str(), ios::in | ios::binary);
  if (!ifs){
    cout << "Error: can't open " << GetDbFilename() << ".bas." <<endl;
    exit(1);
  }
  int temp_i = 0;
  ifs.read(reinterpret_cast<char*>(&temp_i), sizeof(int));
  SetHashSize(temp_i);
  ifs.read(reinterpret_cast<char*>(&temp_i), sizeof(int));
  SetRepeatFlag(temp_i);
  ifs.read(reinterpret_cast<char*>(&temp_i), sizeof(int));
  SetMaximalSpan(temp_i);
  ifs.read(reinterpret_cast<char*>(&temp_i), sizeof(int));
  SetMinAccessibleLength(temp_i);
  ifs.close();
}
