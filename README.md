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
                
    Options:
    (Required)
        -i STR    an RNA sequence in single FASTA format
        -o STR    Output file name
        -d STR    The database name
        
    (Optional)
        -l INT    Max size of seed length [default:20]
        -e INT    Interactin energy threshold for seed search [default: -6.5]
        -f INT    Hybridization energy threshold for removal of the interaction candidate before gapped extension [default: -3.05]
        -x INT    DropOut Length in gapped extension [defualt:18]
        -y INT    DropOut Length in gappless extension [defualt:5]

##Output file format
    Query:ENST00000448587.1|ENSG00000223573.2|OTTHUMG00000150390.2|OTTHUMT00000317918.1|TINCR-001|TINCR|3733|
    Target:ENST00000263461.6|ENSG00000120008.11|OTTHUMG00000019171.4|OTTHUMT00000050707.2|WDR11-001|WDR11|4732|UTR5:1-246|CDS:247-3921|UTR3:3922-4732|
    Energy:-7.41945
    (1274,4730) (1275,4729) (1276,4728) (1277,4727) (1278,4726) (1279,4725) (1280,4724) (1281,4723) 

    Query:ENST00000448587.1|ENSG00000223573.2|OTTHUMG00000150390.2|OTTHUMT00000317918.1|TINCR-001|TINCR|3733|
    Target:ENST00000263461.6|ENSG00000120008.11|OTTHUMG00000019171.4|OTTHUMT00000050707.2|WDR11-001|WDR11|4732|UTR5:1-246|CDS:247-3921|UTR3:3922-4732|
    Energy:-5.9724
    (1344,4730) (1345,4729) (1346,4728) (1347,4727) (1348,4726) (1349,4725) (1350,4724) 

##License
This software is released under the MIT License, see LICENSE.txt.

##Changelogs
2016/08/31 Version 1.0.0 was released.

## Reference
Tsukasa Fukunaga and Michiaki Hamada. "RIblast: An ultrafast RNA-RNA interaction prediction system for comprehensive lncRNA interaction analysis" (in prepearation)
