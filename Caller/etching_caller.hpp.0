//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#ifndef ETCHING_CALLER_HPP
#define ETCHING_CALLER_HPP

#include <set>
#include <map>
#include <unordered_map>
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

const std::string version ;
const std::string PROGRAM = "Program: etching_caller";
const std::string VERSION = "Version: " + version;
const std::string CONTACT = "Contact:\n\tJang-il Sohn (sohnjangil@gmail.com)\n\tJin-Wu Nam (jwnam@hanyang.ac.kr)";

using Position = std::pair < int , int > ;

void print_version();
void caller_usage();

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
		 const bool scanall);
  
void bp_to_vcf ( const std::string genome,
		 const std::string input_bam,
		 const std::string bp_file,
		 const std::string prefix, 
		 const int insert_size, 
		 const bool typing, 
		 const std::string read_orientation, 
		 const double tumor_purity, 
		 const double sequencing_coverage, 
		 const std::string data_type);

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
