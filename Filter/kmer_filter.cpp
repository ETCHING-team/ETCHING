//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#ifndef KMER_FILTER
#define KMER_FILTER

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include <unistd.h>

#include "CPU_TIME_measure.hpp"

void filter_usage(){
  std::cout << "Usage: kmer_filter [options]\n"
	    << "\n"
	    << "Required:\n"
	    << "\t" << "-i <string>\t" << "List of sample (tumor) fastq (required; gz supported).\n"
	    << "\n"
	    << "Filters:\n"
	    << "\t" << "-l <string>\t" << "List of control (matched normal) fastq as (normal) filters (gz supported).\n"
	    << "\t" << "-r <string>\t" << "List of reference genome files to be used as (reference) filters (gz supported).\n"
	    << "\t" << "-a <string>\t" << "Existing filter to add.\n"
	    << "\t" << "           \t" << "If you input PGK, then there must be PGK.kmc_pre and PGK.kmc_suf.\n"
	    << "\n"
	    << "Optional:\n"
	    << "\t" << "-o <string>\t" << "Output file name [filtered_kmer]\n"
	    << "\t" << "-c <int>   \t" << "k-mer histogram cut-off for remove sequencing error [automatically calculated].\n"
	    << "\t" << "-T <string>\t" << "Sequencing data type. W for WGS, or P for Panel [W]\n"
	    << "\t" << "           \t" << "If PANEL, must be specified.\n"
	    << "\t" << "-m <int>   \t" << "Maximum RAM memory in GB [12].\n"
	    << "\t" << "-t <int>   \t" << "Number of thread [8].\n"
	    << "\n"
	    << "Contact:\n\tJang-il Sohn (sohnjangil@gmail.com)\n\tJin-Wu Nam (jwnam@hanyang.ac.kr)\n";
}

int main ( int argc , char ** argv ){
  if ( argc == 1 ){
    filter_usage();
    return 0;
  }

  // const double walltime0(get_wall_time());

  std::string control_list_file;
  std::string reference_list_file;
  std::string outfile="filtered_kmer";

  int control_count = 0;
  int reference_count = 0;

  std::string infile_list;
  std::vector < std::string > file_list;
  std::string existing;

  int maxRAM = 12;
  int num_threads = 8 ;

  int cutoff = 0 ;
  int kl = 31;

  std::string data_type="W";

  int opt;
  while ((opt = getopt ( argc , argv , "i:o:l:r:a:m:t:T:c:k:h") ) != -1 ){
    switch(opt){
    case 'i': infile_list = optarg; break;
    case 'o': outfile = optarg; break;
    case 'l': control_list_file = optarg; break;
    case 'r': reference_list_file = optarg; break;
    case 'a': existing = optarg; break;
    case 'T': data_type = optarg; break;
    case 'm': maxRAM = atoi ( optarg ); break;
    case 't': num_threads = atoi ( optarg ); break;
    case 'c': cutoff = atoi ( optarg ); break;
    case 'k': kl = atoi (optarg ); break;
    case 'h': filter_usage() ; return 0;
    default: std::cerr << "ERROR!!! Invalid option\n" ; return 0;
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////
  //
  // Option check
  //

  if ( data_type != "W" && data_type != "P"){
    std::cerr << "ERROR!!! -T must be one of W, or P. Default is W.\n";
    std::cerr << "\t" << "-T <string>\t" << "Sequencing data type. W for WGS, or P for Panel [W]" ;
    return 0;
  }

  if ( infile_list.size() == 0 ){
    std::cerr << "ERROR!!! -i is required.\n";
    return 0;
  }
  else {
    std::string tmp;
    std::ifstream fin ( infile_list );
    while ( fin >> tmp ){
      if ( tmp.substr(tmp.size()-6) != ".fastq" && 
	   tmp.substr(tmp.size()-3) != ".fq" && 
	   tmp.substr(tmp.size()-9) != ".fastq.gz" && 
	   tmp.substr(tmp.size()-6) != ".fq.gz"){
	std::cerr << "ERROR: Please check files in " << tmp <<"\n";
	std::cerr << "This file must be fastq format.\n";
      }
      file_list.push_back(tmp);
    }
    fin.close();
  }
  
  if ( control_list_file.size() != 0 ){
    std::string tmp;
    std::ifstream fin ( control_list_file.c_str() );
    while ( fin >> tmp){
      if ( tmp.substr(tmp.size()-6) != ".fastq" && 
	   tmp.substr(tmp.size()-3) != ".fq" && 
	   tmp.substr(tmp.size()-9) != ".fastq.gz" && 
	   tmp.substr(tmp.size()-6) != ".fq.gz"){
	std::cerr << "ERROR: Please check files in " << tmp <<"\n";
	std::cerr << "This file must be fastq format.\n";
      }
      file_list.push_back(tmp); 
      control_count ++ ;
    }
    fin.close();
  }
  

  if ( reference_list_file.size() != 0 ){
    std::string tmp;
    std::ifstream fin ( reference_list_file.c_str() );
    while ( fin >> tmp){
      if ( tmp.substr(tmp.size()-6) != ".fasta" && 
	   tmp.substr(tmp.size()-3) != ".fa" && 
	   tmp.substr(tmp.size()-9) != ".fasta.gz" && 
	   tmp.substr(tmp.size()-6) != ".fa.gz" && 
	   tmp.substr(tmp.size()-7) != ".fna.gz" && 
	   tmp.substr(tmp.size()-4) != ".fna" ){
	std::cerr << "ERROR: Please check files in " << reference_list_file <<"\n";
	std::cerr << "The files in this list must be fasta format.\n";
	return 0;
      }
      file_list.push_back(tmp);
      reference_count ++ ;
    }
    fin.close();
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  //
  // File check
  //

  for ( auto & i : file_list ){
    std::ifstream fin (i);
    if ( ! fin ){
      std::cerr << "ERROR!!! There is no " + i + ".\n";
      return 0;
    }
    fin.close();
  }

  if ( existing.size() != 0 ){
    std::string tmp = existing + ".kmc_pre";
    std::ifstream fin;
    fin.open( tmp.c_str() );
    if ( ! fin ){
      std::cerr << "ERROR!!! There is no " + existing+ ".kmc_pre.\n";
      return 0;
    }
    fin.close();
    tmp = existing + ".kmc_suf";
    fin.open( tmp.c_str() );
    if ( ! fin ){
      std::cerr << "ERROR!!! There is no " + existing+ ".kmc_pre.\n";
      return 0;
    }
  }

  if ( file_list.size() + existing.size() == 0){
    std::cerr << "ERROR!!! You did not input any filters.\n";
    std::cerr << "At least one of -l, -r, or -a is required.\n";
    return 0;
  }

  ///////////////////////////////////////////////////////////////////////////////////////
  //
  // Making tmp dir
  //
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::mt19937 generator (seed);

  std::string tmp_dir = "tmp_" + std::to_string(generator());
  std::string command = "mkdir " + tmp_dir;
  std::cout << command << "\n" ; 
  system(command.c_str());

  ///////////////////////////////////////////////////////////////////////////////////////
  //
  // Making normal filter
  //
  //std::cout << "[Making filter]\n";
  if ( control_count !=0 ){
    //std::cout << "Making control filter\n";
    command = "kmc -v -k" + std::to_string(kl) + " -t"+ std::to_string ( num_threads ) + " -m" + std::to_string ( maxRAM ) + " -ci2 -fq @" + control_list_file + " control_filter " + tmp_dir + " > control_filter.log 2> control_filter.err";
    std::cout << command << "\n" ; 
    system(command.c_str());
  }
  ///////////////////////////////////////////////////////////////////////////////////////
  //
  // Making reference filter
  //
  if ( reference_count != 0){
    //std::cout << "Making reference filter\n";
    command = "kmc -v -k" + std::to_string(kl) + " -t"+ std::to_string ( num_threads ) + " -m" + std::to_string ( maxRAM ) + " -ci1 -fm @" + reference_list_file + " reference_filter " + tmp_dir + " > reference_filter.log 2> reference_filter.err";
    std::cout << command << "\n" ; 
    system(command.c_str());
  }
  ///////////////////////////////////////////////////////////////////////////////////////
  //
  // Merging filters
  //
  if ( existing.size() != 0 ){
    //std::cout << "Merging filters\n"; 
    if ( control_count !=0 && reference_count == 0 ){
      command = "kmc_tools simple control_filter " + existing + " union filter 2> merging.err";
      std::cout << command << "\n";
      system(command.c_str() );
    }
    else if ( control_count ==0 && reference_count != 0 ){
      command = "kmc_tools simple reference_filter " + existing + " union filter";
      std::cout << command << "\n" ;
      system(command.c_str() );
    }
    else if ( control_count !=0 && reference_count != 0 ){
      command = "kmc_tools simple control_filter " + existing + " union filter_tmp 2> merging_tmp.err";
      std::cout << command << "\n";
      system(command.c_str() );
      command = "kmc_tools simple filter_tmp reference_filter union filter 2> merging.err";
      std::cout << command << "\n";
      system(command.c_str() );
    }
    else {
      command = "ln -sf " + existing + ".kmc_pre filter.kmc_pre";
      std::cout << command << "\n"; 
      system(command.c_str() );
      command = "ln -sf " + existing + ".kmc_suf filter.kmc_suf";
      std::cout << command << "\n"; 
      system(command.c_str() );
    }
  }
  else{
    if ( control_count !=0 && reference_count == 0 ){
      command = "ln -sf control_filter.kmc_pre filter.kmc_pre";
      std::cout << command << "\n"; 
      system(command.c_str() );
      command = "ln -sf control_filter.kmc_suf filter.kmc_suf";
      std::cout << command << "\n"; 
      system(command.c_str() );
    }
    else if ( control_count ==0 && reference_count != 0 ){
      command = "ln -sf reference_filter.kmc_pre filter.kmc_pre";
      std::cout << command << "\n"; 
      system(command.c_str() );
      command = "ln -sf reference_filter.kmc_suf filter.kmc_suf";
      std::cout << command << "\n"; 
      system(command.c_str() );
    }
    else if ( control_count !=0 && reference_count != 0 ){
      std::cout << "Merging filters\n"; 
      command = "kmc_tools simple control_filter reference_filter union filter > merging.log 2> merging.err";
      std::cout << command << "\n"; 
      system(command.c_str() );
    }
    else {
      std::cerr <<"ERROR!!! No filter was input";
      return 0;
    }
  }
  // const double walltime1(get_wall_time());
  // std::cout << "Wall-clock time: " << walltime1 - walltime0 << " sec" << "\n";
  // std::cout << "\n";

  ///////////////////////////////////////////////////////////////////////////////////////
  //
  // Making tumor k-mer table
  //
  //std::cout << "[ETCHING filter]\n";
  //std::cout << "Making Sample k-mer DB\n";
  command = "kmc -v -k" + std::to_string(kl) + " -t"+ std::to_string ( num_threads ) + " -m" + std::to_string ( maxRAM ) + " -ci2 -fq @" + infile_list + " sample " + tmp_dir + " > sample_kmer_table.log 2> sample_kmer_table.err";
  std::cout << command << "\n" ; 
  system(command.c_str());

  ///////////////////////////////////////////////////////////////////////////////////////
  //
  // Calculating cut-off
  //
  if ( cutoff <= 0 ){
    std::vector <std::size_t> depth;
    std::vector <std::size_t> freq;
    std::size_t d,f;
    std::size_t local_min;
  
    std::string sample_hist = "sample.hist";

    command = "kmc_tools transform sample histogram " + sample_hist + " 2> sample_kmer_hist.err";
    std::cout << command << "\n"; 
    system(command.c_str() );

    std::ifstream fin(sample_hist.c_str());
    while ( fin >> d >> f ){
      depth.push_back(d);
      freq.push_back(f);
    }
    fin.close();

    if ( data_type == "W" ){ // for WGS
      for ( size_t i = 1 ; i < depth.size() ; i ++ ){
	if ( freq[i-1] < freq[i] ){
	  local_min = depth[i-1];
	  break;
	}
      }
    }
    else if ( data_type == "P" ) { // for panel
      for ( size_t i = 1 ; i < depth.size() ; i ++ ){
	if ( freq[0] > 10 * freq[i] ){
	  // local_min = depth[i];
	  if ( freq[i-1] - freq[0]/10 < freq[0]/10 - freq[i] ){
	    local_min = depth[i-1];
	  }
	  else{
	    local_min = depth[i];
	  }
	  break;
	}
      }
    }
    if ( local_min < 2 ) local_min = 2 ;
    cutoff = local_min;
  }
  //std::cout << "k-mer cutoff: " << cutoff << "\n";


  ///////////////////////////////////////////////////////////////////////////////////////
  //
  // Applying filter
  //
  //std::cout << "Appying filter\n";
  command = "kmc_tools simple sample -ci" + std::to_string(cutoff) + " filter -ci1 kmers_subtract " + outfile + " 2> " + outfile + ".err";
  std::cout << command << "\n" ;
  system(command.c_str() );

  ///////////////////////////////////////////////////////////////////////////////////////
  //
  // Print filter
  //
  //std::cout << "Writing sample-specific k-mers\n";
  command = "kmc_dump " + outfile + " " + outfile ;
  std::cout << command << "\n"; 
  system(command.c_str() );


  ///////////////////////////////////////////////////////////////////////////////////////
  //
  // End
  //

  command = "rm -rf " + tmp_dir ;
  std::cout << command << "\n" ;
  system(command.c_str());

  // const double walltime2(get_wall_time());
  //std::cout << "[Finished]\n";
  //std::cout << "Wall-clock time: " << walltime2 - walltime0 << " sec" << "\n";

  return 0;
}


#endif
