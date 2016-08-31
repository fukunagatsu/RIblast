# RIblast
RIblast is ultrafast RNA-RNA interaction prediction software for comprehensive lncRNA interactome analysis based on seed-and-extension algorithm

##Version
Version 1.0.0 (2016/08/31)

##Acknowledgements
We used BL* energy model as RNA secondary structure energy model.
You can download this energy model from  http://www.cs.ubc.ca/labs/beta/Projects/RNA-Params
We thank Dr. Mirela Andronescu for the development of this model.

##Usage
RIblast consists of database construction search step and rna interactin search step. Firstly, you have to generate formatted database files from FASTA fomatted RNA sequences file using RIblast "db" command. RIblast "db" command requires 2 options ([-i InputFastaFile] and [-o OutputDbName]). Then, you can search RNA-RNA interaction of a query seuence to database sequences using RIblast "ris" command. RIblast "ris" command requires 3 options ([-i InputFastaFile], [-o OutputFileName] and [-d DatabaseFileName]).

##Example
    ./RIblast db -i dbRNA.fa -o test_db
    ./RIblast ris -i queryRNA.fa -o output.txt -d test_db

##Command and Options
    db: convert a FASTA file to RIblast database files  

     RIblast db [-i InputFastaFile] [-o OutputDbName] [-r RepeatMaskingStyle]  
                [-s LookupTableSize] [-w MaximalSpan] [-d MinAccessibleLength]  
    Options:
    (Required)
        -i STR    RNA sequences in FASTA format
        -o STR    The database name
        
    (Optional) 
        -r INT    Designation of repeat masking style 0:hard-masking, 1:soft-masking, 2:hard-masking [default:0]
        -s INT    Lookup table size of short string search [default: 8]
        -w INT    The constraint of maximal distance between the bases that form base pairs [default: 70]
        -d INT    Minimum accessible length for accessibility approximation [defualt:5]
        
    ris: search RNA-RNA interaction between a query and database sequences
    
    RIblast ris [-i InputFastaFile] [-o OutputFileName] [-d DatabaseFileName]
                [-l MaxSeedLength] [-e HybridizationEnergyThreshold] [-f InteractionEnergyThreshold]
                [-x DropOutLengthInGappedExtension] [-y DropOutLengthInUngappedExtension]
                
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

##License
This software is released under the MIT License, see LICENSE.txt.

##Changelogs
2016/08/31 Version 1.0.0 was released.

## Reference
Tsukasa Fukunaga and Michiaki Hamada. "RIblast: An ultrafast RNA-RNA interaction prediction system for comprehensive lncRNA interaction analysis" (in prepearation)
