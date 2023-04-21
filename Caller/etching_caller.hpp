//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#ifndef ETCHING_CALLER_HPP
#define ETCHING_CALLER_HPP

#include <set>
#include <map>
#include <unordered_set>
#include <iostream>
#include <fstream>
#include <iostream>
#include <string>
#include <string.h>
#include <unistd.h>
#include <limits>
#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>

#include <api/BamAlignment.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>

#include "my_vcf.hpp"

#include "Peak_memory_measure.hpp"
#include "CPU_TIME_measure.hpp"

using Position = std::pair < int , int > ;

void caller_usage(int argc , char ** argv);

template < class T >
class Format{
private:
  std::map < std::string , T > key;
public:
  Format(){
  }
  ~Format(){};
  
  T & operator [](std::string cutoff){
    return key[cutoff];
  }
};


void find_path ( const std::string input_bam, 
		 const std::string prefix, 
		 const int insert_size, 
		 const bool scanall,
		 const bool rescue );
  
void bp_to_vcf ( const std::string genome,
		 const std::string input_bam,
		 const std::string bp_file,
		 const std::string prefix, 
		 const int insert_size, 
		 const bool typing, 
		 const std::string read_orientation, 
		 const double tumor_purity, 
		 const double sequencing_coverage);


int return_discordant_pair_number ( BamTools::BamReader & reader1, 
				    const Position Pos1, 
				    const Position Pos2, 
				    const int insert_size, 
				    BamTools::RefVector & references);

int return_discordant_pair_number ( BamTools::BamReader & reader1, 
				    const Position Pos1, 
				    const int insert_size, 
				    BamTools::RefVector & references);

double get_wall_time();
double get_cpu_time();


std::size_t get_peak_memory();

#endif
