#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include "db_construction.h"
#include "db_construction_parameters.h"
#include "rna_interaction_search.h"
#include "rna_interaction_search_parameters.h"

using namespace std;

void PrintUsage() {
  cout << "RIblast - RNA-RNA interaction preriction tool. version 1.0" << endl;
  cout << "\n";
  cout << "Options\n";
  cout << "db: convert a FASTA file to RIblast database files\n";
  cout << "\n";
  cout << "RIblast db [-i InputFastaFile] [-o OutputDbName] [-r RepeatMaskingStyle]\n";
  cout << "           [-s LookupTableSize] [-w MaximalSpan] [-d MinAccessibleLength]\n";
  cout << "\n";
  cout << "  Options:\n";
  cout << " (Required)\n";
  cout << "    -i STR    RNA sequences in FASTA format\n";
  cout << "    -o STR    The database name\n";
  cout << "\n";
  cout << " (Optional)\n";
  cout << "    -r INT    Designation of repeat masking style 0:hard-masking, 1:soft-masking, 2:hard-masking [default:0]\n";
  cout << "    -s INT    Lookup table size of short string search [default: 8]\n";
  cout << "    -w INT    The constraint of maximal distance between the bases that form base pairs [default: 70]\n";
  cout << "    -d INT    Minimum accessible length for accessibility approximation [defualt:5]\n";
  cout << "\n";
  cout << "\n";
  cout << "ris: search RNA-RNA interaction between a query and database sequences\n";
  cout << "\n";
  cout << "RIblast ris [-i InputFastaFile] [-o OutputFileName] [-d DatabaseFileName]\n";
  cout << "            [-l MaxSeedLength] [-e HybridizationEnergyThreshold] [-f InteractionEnergyThreshold]\n";
  cout << "            [-x DropOutLengthInGappedExtension] [-y DropOutLengthInUngappedExtension]\n";
  cout << "\n";
  cout << "  Options:\n";
  cout << " (Required)\n";
  cout << "    -i STR    an RNA sequence in single FASTA format\n";
  cout << "    -o STR    Output file name\n";
  cout << "    -d STR    The database name\n";
  cout << "\n";
  cout << " (Optional)\n";
  cout << "    -l INT    Max size of seed length [default:20]\n";
  cout << "    -e INT    Interactin energy threshold for seed search [default: -6.5]\n";
  cout << "    -f INT    Hybridization energy threshold for removal of the interaction candidate before gapped extension [default: -3.05]\n";
  cout << "    -x INT    DropOut Length in gapped extension [defualt:18]\n";
  cout << "    -y INT    DropOut Length in gappless extension [defualt:5]\n";
}

int main(int argc, char* argv[]) {
  if(argc == 1 || strcmp(argv[1],"h") == 0){
    PrintUsage();
    exit(1);
  }
  
  if(strcmp(argv[1],"db") == 0){
    DbConstructionParameters parameters;
    parameters.SetParameters(argc - 1, argv + 1);
    DbConstruction db_construction;
    db_construction.Run(parameters);
    
  }else if(strcmp(argv[1],"ris") == 0){
    RnaInteractionSearchParameters parameters;
    parameters.SetParameters(argc - 1, argv + 1);
    parameters.SetDbParameters();
    RnaInteractionSearch rna_interaction_search(parameters.GetHashSize());
    rna_interaction_search.Run(parameters);
  }else{
    cerr << "Error: You must specify the mode of RIblast (db or ris)." << endl;
    exit(1);
  }
  return(0);
}
