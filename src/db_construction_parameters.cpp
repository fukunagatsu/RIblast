/*
 *  db_construction_parameters.cpp
 *
 *  Created on: 2016/8/31
 *      Author: Tsukasa Fukunaga
 */

#include "db_construction_parameters.h"
#include <getopt.h>
#include <stdlib.h>

void DbConstructionParameters::SetParameters(int argc, char* argv[]) {
  int c;
  extern char *optarg;
  while ((c = getopt(argc, argv, "i:o:r:s:w:d:t:")) != -1) {
    switch (c) {
    case 'i':
      _input_filename = optarg;
      break;

    case 'o':
      _db_filename = optarg;
      break;

    case 'r':
      _repeat_flag = atoi(optarg);
      break;

    case 's':
      _hash_size = atoi(optarg);
      break;

    case 'w':
      _maximal_span = atoi(optarg);
      break;
      
    case 'd':
      _min_accessible_length = atoi(optarg);
      break;

    default:
      cerr << "The argument is invalid command." << endl;
      exit(1);
    }
  }
}
