/*
 * fastafile_reader.cpp
 *
 *  Created on: 2016/8/31
 *      Author: Tsukasa Fukunaga
 */

#include "fastafile_reader.h"
#include <fstream>
#include <iostream>

void FastafileReader::ReadFastafile(string input_file_name, vector<string> &sequences, vector<string> &names){
  ifstream fp;
  string buffer;
  fp.open(input_file_name.c_str(),ios::in);
  if (!fp){
    cout << "Error: can't open input_file:"+input_file_name+"." <<endl;
    exit(1);
  }
  getline(fp,buffer);
  names.push_back(buffer.substr(1,buffer.size()-1));
  string temp_sequence = "";
  int count = 0;
  while(getline(fp,buffer)){
    if(buffer[0] == '>'){
      names.push_back(buffer.substr(1,buffer.size()-1));
      sequences.push_back(temp_sequence);
      temp_sequence = "";
    }else{
      if(buffer.size()>=2){
	if(buffer.substr(buffer.size()-2,2) == "\r\n"){
	  buffer.erase(buffer.size()-2,2);
	}
      }
      if(buffer[buffer.size()-1] == '\r' || buffer[buffer.size()-1] == '\n'){
	buffer.erase(buffer.size()-1,1);
      }
      count += 1;
      temp_sequence = temp_sequence + buffer;
    }
  }
  sequences.push_back(temp_sequence);
  fp.close();
}

void FastafileReader::ReadFastafile(string input_file_name, string &sequence, string &name){
  ifstream fp;
  string buffer;
  fp.open(input_file_name.c_str(),ios::in);
  if (!fp){
    cout << "Error: can't open input_file:"+input_file_name+"." <<endl;
    exit(1);
  }
  getline(fp,buffer);
  name = buffer.substr(1,buffer.size()-1);
  sequence = "";
  while(getline(fp,buffer)){
    if(buffer.substr(buffer.size()-2,2) == "\r\n"){
      buffer.erase(buffer.size()-2,2);
    }else if(buffer[buffer.size()-1] == '\r' || buffer[buffer.size()-1] == '\n'){
      buffer.erase(buffer.size()-1,1);
    }
    sequence = sequence + buffer;
  }
  fp.close();
}
