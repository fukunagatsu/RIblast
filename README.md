# RIblast
RIblast is ultrafast RNA-RNA interaction prediction software based on seed-and-extension algorithm for comprehensive lncRNA interactome analysis.

##Version
Version 1.0.0 (2016/08/31)

##Acknowledgements
We used BL* energy model as RNA secondary structure energy model.
You can download this energy model from  http://www.cs.ubc.ca/labs/beta/Projects/RNA-Params
We thank Dr. Mirela Andronescu for the development of this model.

##Usage
RIblast consists of database construction search step and rna interactin search step. Firstly, you have to generate formatted database files from FASTA fomatted RNA sequences file using RIblast "db" command. RIblast "db" command requires 2 options ([-i InputFastaFile] and [-o OutputDbName]). Then, you can search RNA-RNA interaction of a query seuence to database sequences using RIblast "ris" command. RIblast "ris" command requires 3 options ([-i InputFastaFile], [-o OutputFileName] and [-d DatabaseFileName]). The number of sequence in query fasta file must be 1.

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
        -r INT    Designation of repeat masking style 0:hard-masking, 1:soft-masking, 2:no-masking [default:0]
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
        -e INT    Hybridization energy threshold for seed search [default: -6.5]
        -f INT    Interaction energy threshold for removal of the interaction candidate before gapped extension [default: -3.05]
        -x INT    DropOut Length in gapped extension [defualt:18]
        -y INT    DropOut Length in ungapped extension [defualt:5]

##Output file format
RIblast outputs detected interactions as follows.
An interaction is expressed in four lines, and blank line is inserted between each interaction.
The first, second, third and fourth line of an interaction describes name of query RNA, name of target RNA, interaction energy of the interaction, and interacted positions ([position in query],[position in target]), respectively.

Example  
Query:qrna  
Target:taget_rna1  
Energy:-7.41945  
(4,30) (5,29) (6,28) (7,27) (8,26) (9,25) (10,24) (11,23)  

Query:qrna  
Target:taget_rna2  
Energy:-9.9724  
(72,185) (73,184) (74,183) (77,176) (78,175) (79,174) (80,173) (81,172) (82,171) (83,170) 

##License
This software is released under the MIT License, see LICENSE.txt.

##Changelogs
2016/08/31 Version 1.0.0 was released.

## Reference
Tsukasa Fukunaga and Michiaki Hamada. "RIblast: An ultrafast RNA-RNA interaction prediction system for comprehensive lncRNA interaction analysis" (in prepearation)
