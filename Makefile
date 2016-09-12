CXXFLAGS = -O3

RIblast: main.cpp db_construction_parameters.cpp db_construction.cpp  fastafile_reader.cpp encoder.cpp sais.c raccess.cpp energy_par.h rna_interaction_search.cpp rna_interaction_search_parameters.cpp seed_search.cpp ungapped_extension.cpp gapped_extension.cpp

	  $(CXX) $(CXXFLAGS) -o RIblast main.cpp db_construction_parameters.cpp fastafile_reader.cpp db_construction.cpp encoder.cpp raccess.cpp rna_interaction_search.cpp rna_interaction_search_parameters.cpp ungapped_extension.cpp seed_search.cpp gapped_extension.cpp -x c sais.c
