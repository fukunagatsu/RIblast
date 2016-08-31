#ifndef FASTAFILE_READER_H
#define FASTAFILE_READER_H

#include <string>
#include <vector>
#include <stdlib.h>

using namespace std;

class FastafileReader {
 public:
  FastafileReader() {}
  void ReadFastafile(string input_file_name, vector<string> &sequences, vector<string> &names);
  void ReadFastafile(string input_file_name, string &sequences, string &name);
};

#endif
