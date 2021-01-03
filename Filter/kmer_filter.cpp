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
	    << "\t" << "-i <string>\t" << "List of sample (tumor) in fastq(.gz) or bam.\n"
	    << "\t" << "           \t" << "If the file format is fastq, you must have a pair of fastq files of paired-end reads.\n"
	    << "\t" << "           \t" << "If the file format is bam, use only one mapped bam file of paired-end reads\n"
	    << "\t" << "           \t" << "Do not use fastq and bam files at the same time.\n"
	    << "\n"
	    << "Filters:\n"
	    << "\t" << "-c <string>\t" << "List of control (matched normal) in fastq(.gz) or bam (mapped or unmapped).\n"
	    << "\t" << "-r <string>\t" << "List of control (reference genome) files in fasta(.gz).\n"
	    << "\t" << "-a <string>\t" << "Existing filter to add. If you input PGK, you must have PGK.kmc_pre and PGK.kmc_suf.\n"
	    << "\n"
	    << "Options:\n"
	    << "\t" << "-p <string>\t" << "Output file name [filtered_kmer]\n"
	    << "\t" << "-l <string>\t" << "k-mer size [31]\n"
	    << "\t" << "-K <int>   \t" << "k-mer histogram cut-off for remove sequencing error [automatically calculated].\n"
	    << "\t" << "-M <int>   \t" << "Maximum k-mer frequency to be counted [10000].\n"
	    << "\t" << "-T <string>\t" << "Sequencing data type. W for WGS, or P for Panel [W]\n"
	    << "\t" << "           \t" << "If PANEL, must be specified.\n"
	    << "\t" << "-m <int>   \t" << "Maximum RAM memory in GB [12].\n"
	    << "\t" << "-t <int>   \t" << "Number of thread [8].\n"
	    << "\t" << "-D <int>   \t" << "Directory of KMC-3 [null].\n"
	    << "\t" << "-E         \t" << "Store KMC database files"
	    << "\n"
	    << "Contact:\n\tJang-il Sohn (sohnjangil@gmail.com)\n\tJin-Wu Nam (jwnam@hanyang.ac.kr)\n";
}

void print_file ( std::string infile ){
  std::ifstream fin(infile.c_str());
  std::string tmp;
  std::cout << infile << "\n";
  while ( fin >> tmp ){
    std::cout << "\t" << tmp << "\n";
  }
}


int main ( int argc , char ** argv ){
  if ( argc == 1 ){
    filter_usage();
    return 0;
  }

  // const double walltime0(get_wall_time());

  std::string control_list_file;
  std::string reference_list_file;
  std::string prefix="filtered_kmer";

  int64_t control_count = 0;
  int64_t reference_count = 0;

  std::string infile_list;
  std::vector < std::string > file_list;
  std::string existing;

  int64_t maxRAM = 12;
  int64_t num_threads = 8 ;

  int64_t cutoff = 0 ;
  int64_t maxkfreq = 10000;
  int64_t kl = 31;

  std::string data_type="W";

  std::string DIR;

  bool store_kmc = 0 ;

  int64_t opt = 0;
  while ((opt = getopt ( argc , argv , "i:p:c:r:a:m:t:T:K:M:l:D:Eh") ) != -1 ){
    switch(opt){
    case 'i': infile_list = optarg; break;
    case 'p': prefix = optarg; break;
    case 'c': control_list_file = optarg; break;
    case 'r': reference_list_file = optarg; break;
    case 'a': existing = optarg; break;
    case 'T': data_type = optarg; break;
    case 'm': maxRAM = atoi ( optarg ); break;
    case 't': num_threads = atoi ( optarg ); break;
    case 'K': cutoff = atoi ( optarg ); break;
    case 'M': maxkfreq = atoi ( optarg ); break;
    case 'l': kl = atoi (optarg ); break;
    case 'D': DIR = optarg ; break;
    case 'E': store_kmc = 1; break;
    case 'h': filter_usage() ; return 0;
    default: std::cout << "ERROR!!! Invalid option: " << optarg << "\n" ; return -1;
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////
  //
  // Option check
  //

  // pop_back if DIR has "/" charactor at the tail
  if ( DIR.size() != 0 ){
    if ( DIR[DIR.size()-1] != '/' ){
      DIR += '/';
    }
  }

  if ( data_type != "W" && data_type != "P"){
    std::cout << "ERROR!!! -T must be one of W, or P. Default is W.\n";
    std::cout << "\t" << "-T <string>\t" << "Sequencing data type. W for WGS, or P for Panel [W]" ;
    return -1;
  }

  if ( infile_list.size() == 0 ){
    std::cout << "ERROR!!! -i is required.\n";
    return -1;
  }
  else {
    std::string tmp;
    std::ifstream fin ( infile_list );

    while ( fin >> tmp ){
      if ( tmp.substr(tmp.size()-6) != ".fastq" && 
	   tmp.substr(tmp.size()-3) != ".fq" && 
	   tmp.substr(tmp.size()-9) != ".fastq.gz" && 
	   tmp.substr(tmp.size()-6) != ".fq.gz" &&
	   tmp.substr(tmp.size()-4) != ".bam" ){
	std::cout << "ERROR!!! Please check files in " << tmp <<".\n";
	std::cout << "This file must be fastq(.gz) or bam format.\n";
      }
      //count ++;
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
	   tmp.substr(tmp.size()-6) != ".fq.gz" &&
	   tmp.substr(tmp.size()-4) != ".bam" ){
	std::cout << "ERROR: Please check files in " << tmp <<".\n";
	std::cout << "This file must be fastq(.gz) or bam format.\n";
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
	std::cout << "ERROR: Please check files in " << reference_list_file <<".\n";
	std::cout << "The files in this list must be fasta(.gz) format.\n";
	return -1;
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
      std::cout << "ERROR!!! There is no " + i + ".\n";
      return -1;
    }
    fin.close();
  }

  if ( existing.size() != 0 ){
    std::string tmp = existing + ".kmc_pre";
    std::ifstream fin;
    fin.open( tmp.c_str() );
    if ( ! fin ){
      std::cout << "ERROR!!! There is no " + existing+ ".kmc_pre.\n";
      return -1;
    }
    fin.close();
    tmp = existing + ".kmc_suf";
    fin.open( tmp.c_str() );
    if ( ! fin ){
      std::cout << "ERROR!!! There is no " + existing+ ".kmc_pre.\n";
      return -1;
    }
  }

  if ( file_list.size() + existing.size() == 0){
    std::cout << "ERROR!!! You did not input any filters.\n";
    std::cout << "At least one of -l, -r, or -a is required.\n";
    return -1;
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

  std::vector < std::string > bam_files ;
  std::vector < std::string > fastq_files ;
  std::string file_format;

  int64_t bam_count = 0;

  if ( control_count !=0 ){

    if ( control_list_file.size() != 0 ){
      std::string tmp;
      std::ifstream fin ( control_list_file.c_str() );
      
      int64_t fastq_count = 0;
      
      while ( fin >> tmp){
	std::cout << "[Filter as control : " << tmp << "]\n";
	if ( tmp.substr(tmp.size()-4) == ".bam" ){
	  bam_count ++;
	  bam_files.push_back(tmp);
	}
	else{
	  fastq_count ++;
	  fastq_files.push_back(tmp);
	}
      }
      fin.close();

      if ( bam_count == 0 && fastq_count != 0 ){
	file_format = "q";
	command 
	  = DIR + "kmc -v -k" + std::to_string(kl) + " -t"+ std::to_string ( num_threads ) + " -m" + std::to_string ( maxRAM ) 
	  + " -ci2 -f" + file_format + " @" + control_list_file + " control_filter " + tmp_dir + " 2>&1";
	std::cout << command << "\n" ; 
	print_file(control_list_file);
	system(command.c_str());
	
      }
      else if ( bam_count != 0 && fastq_count == 0 ){
	file_format = "bam";
	command 
	  = DIR + "kmc -v -k" + std::to_string(kl) + " -t"+ std::to_string ( num_threads ) + " -m" + std::to_string ( maxRAM ) 
	  + " -ci2 -f" + file_format + " @" + control_list_file + " control_filter " + tmp_dir + " 2>&1";
	std::cout << command << "\n" ; 
	print_file(control_list_file);
	system(command.c_str());
	
      }
      else if ( fastq_count != 0 && bam_count != 0){

	std::ofstream fout;
	fout.open("control_filter_fastq_file_list.txt");
	for ( auto i : fastq_files ){
	  fout << i << "\n";
	}
	fout.close();
	fastq_files.clear();

	fout.open("control_filter_bam_file_list.txt");
	for ( auto i : bam_files ){
	  fout << i << "\n";
	}
	fout.close();

	file_format = "q";
	command 
	  = DIR + "kmc -v -k" + std::to_string(kl) + " -t"+ std::to_string ( num_threads ) + " -m" + std::to_string ( maxRAM ) 
	  + " -ci2 -f" + file_format + " @control_filter_fastq_file_list.txt control_filter_fastq " + tmp_dir + " 2>&1";
	std::cout << command << "\n" ; 
	print_file("control_filter_fastq_file_list.txt");
	system(command.c_str());

	
	file_format = "bam";
	command 
	  = DIR + "kmc -v -k" + std::to_string(kl) + " -t"+ std::to_string ( num_threads ) + " -m" + std::to_string ( maxRAM ) 
	  + " -ci2 -f" + file_format + " @control_filter_bam_file_list.txt control_filter_bam " + tmp_dir + " 2>&1";
	std::cout << command << "\n" ; 
	print_file("control_filter_bam_file_list.txt");
	system(command.c_str());
	

	command = DIR + "kmc_tools simple control_filter_fastq control_filter_bam union control_filter 2>&1";
	std::cout << command << "\n" ; 
	system(command.c_str());
	if ( store_kmc == 0 ) {
	  std::cout << "remove(\"control_filter_fastq.kmc_pre\")\n";
	  remove("control_filter_fastq.kmc_pre");
	  std::cout << "remove(\"control_filter_fastq.kmc_suf\")\n";
	  remove("control_filter_fastq.kmc_suf");
	  std::cout << "remove(\"control_filter_bam.kmc_pre\")\n";
	  remove("control_filter_bam.kmc_pre");
	  std::cout << "remove(\"control_filter_bam.kmc_suf\")\n";
	  remove("control_filter_bam.kmc_suf");
	}
      }
    }

    //std::cout << "Making control filter\n";
  }
  ///////////////////////////////////////////////////////////////////////////////////////
  //
  // Making reference filter
  //
  if ( reference_count != 0){
    //std::cout << "Making reference filter\n";
    command 
      = DIR + "kmc -v -k" + std::to_string(kl) + " -t"+ std::to_string ( num_threads ) + " -m" + std::to_string ( maxRAM ) 
      + " -ci1 -fm @" + reference_list_file + " reference_filter " + tmp_dir + " 2>&1";
    std::cout << command << "\n" ; 
    print_file(reference_list_file);
    system(command.c_str());
    
  }
  ///////////////////////////////////////////////////////////////////////////////////////
  //
  // Merging filters
  //
  if ( existing.size() != 0 ){
    //std::cout << "Merging filters\n"; 
    if ( control_count !=0 && reference_count == 0 ){
      command = DIR + "kmc_tools simple control_filter " + existing + " union filter 2>&1";
      std::cout << command << "\n";
      
      system(command.c_str() );
      if ( store_kmc == 0 ) {
	std::cout << "remove(\"control_filter.kmc_pre\")\n";
	remove("control_filter.kmc_pre");
	std::cout << "remove(\"control_filter.kmc_suf\")\n";
	remove("control_filter.kmc_suf");
      }
    }
    else if ( control_count ==0 && reference_count != 0 ){
      command = DIR + "kmc_tools simple reference_filter " + existing + " union filter 2>&1";
      std::cout << command << "\n" ;
      
      system(command.c_str() );
      if ( store_kmc == 0 ) {
	std::cout << "remove(\"reference_filter.kmc_pre\")\n";
	remove("reference_filter.kmc_pre");
	std::cout << "remove(\"reference_filter.kmc_suf\")\n";
	remove("reference_filter.kmc_suf");
      }
    }
    else if ( control_count !=0 && reference_count != 0 ){
      command = DIR + "kmc_tools simple control_filter " + existing + " union filter_tmp 2>&1";
      std::cout << command << "\n";
      
      system(command.c_str() );
      if ( store_kmc == 0 ) {
	std::cout << "remove(\"control_filter.kmc_pre\")\n";
	remove("control_filter.kmc_pre");
	std::cout << "remove(\"control_filter.kmc_suf\")\n";
	remove("control_filter.kmc_suf");
      }

      command = DIR + "kmc_tools simple filter_tmp reference_filter union filter 2>&1";
      std::cout << command << "\n";
      
      system(command.c_str() );
      if ( store_kmc == 0 ) {
	std::cout << "remove(\"filter_tmp.kmc_pre\")\n";
	remove("filter_tmp.kmc_pre");
	std::cout << "remove(\"filter_tmp.kmc_suf\")\n";
	remove("filter_tmp.kmc_suf");
	std::cout << "remove(\"reference_filter.kmc_pre\")\n";
	remove("reference_filter.kmc_pre");
	std::cout << "remove(\"reference_filter.kmc_suf\")\n";
	remove("reference_filter.kmc_suf");
      }
    }
    else {
      // command = "ln -sf " + existing + ".kmc_pre filter.kmc_pre";
      // std::cout << command << "\n"; 
      // system(command.c_str() );
      // command = "ln -sf " + existing + ".kmc_suf filter.kmc_suf";
      // std::cout << command << "\n"; 
      // system(command.c_str() );
      command = "mv -f " + existing + ".kmc_pre filter.kmc_pre";
      std::cout << command << "\n"; 
      system(command.c_str() );
      command = "mv -f " + existing + ".kmc_suf filter.kmc_suf";
      std::cout << command << "\n"; 
      system(command.c_str() );
    }
  }
  else{
    if ( control_count !=0 && reference_count == 0 ){
      // command = "ln -sf control_filter.kmc_pre filter.kmc_pre";
      // std::cout << command << "\n"; 
      // system(command.c_str() );
      // command = "ln -sf control_filter.kmc_suf filter.kmc_suf";
      // std::cout << command << "\n"; 
      // system(command.c_str() );
      command = "mv -f control_filter.kmc_pre filter.kmc_pre";
      std::cout << command << "\n"; 
      system(command.c_str() );
      command = "mv -f control_filter.kmc_suf filter.kmc_suf";
      std::cout << command << "\n"; 
      system(command.c_str() );
    }
    else if ( control_count ==0 && reference_count != 0 ){
      // command = "ln -sf reference_filter.kmc_pre filter.kmc_pre";
      // std::cout << command << "\n"; 
      // system(command.c_str() );
      // command = "ln -sf reference_filter.kmc_suf filter.kmc_suf";
      // std::cout << command << "\n"; 
      // system(command.c_str() );
      command = "mv -f reference_filter.kmc_pre filter.kmc_pre";
      std::cout << command << "\n"; 
      system(command.c_str() );
      command = "mv -f reference_filter.kmc_suf filter.kmc_suf";
      std::cout << command << "\n"; 
      system(command.c_str() );
    }
    else if ( control_count !=0 && reference_count != 0 ){
      std::cout << "Merging filters\n"; 
      command = DIR + "kmc_tools simple control_filter reference_filter union filter 2>&1";
      std::cout << command << "\n"; 
      
      system(command.c_str() );
      if ( store_kmc == 0 ) {
	std::cout << "remove(\"control_filter.kmc_pre\")\n";
	remove("control_filter.kmc_pre");
	std::cout << "remove(\"control_filter.kmc_suf\")\n";
	remove("control_filter.kmc_suf");
	std::cout << "remove(\"reference_filter.kmc_pre\")\n";
	remove("reference_filter.kmc_pre");
	std::cout << "remove(\"reference_filter.kmc_suf\")\n";
	remove("reference_filter.kmc_suf");
      }
    }
    else {
      std::cout <<"ERROR!!! No filter was input";
      return -1;
    }
  }
  // const double walltime1(get_wall_time());
  // std::cout << "Wall-clock time: " << walltime1 - walltime0 << " sec" << "\n";
  // 

  ///////////////////////////////////////////////////////////////////////////////////////
  //
  // Making tumor k-mer table
  //
  //std::cout << "[ETCHING filter]\n";
  //std::cout << "Making Sample k-mer DB\n";
  
  file_format.clear();
  
  std::string tmp;
  std::ifstream fin ( infile_list );
  
  int64_t fastq_count = 0 ;
  bam_count = 0 ;
  bam_files.clear();
  while ( fin >> tmp){
    if ( tmp.substr(tmp.size()-4) == ".bam" ){
      bam_count ++;
      bam_files.push_back(tmp);
    }
    else{
      fastq_count ++;
      fastq_files.push_back(tmp);
    }
  }
  fin.close();


  if ( bam_count == 0 && fastq_count != 0 ){
    file_format = "q";
    command 
      = DIR + "kmc -v -k" + std::to_string(kl) + " -t"+ std::to_string ( num_threads ) + " -m" + std::to_string ( maxRAM ) 
      + " -ci2 " + "-cx" + std::to_string(maxkfreq) + " -f" + file_format + " @" + infile_list + " sample " + tmp_dir + " 2>&1";
    std::cout << command << "\n" ; 
    print_file ( infile_list );
    system(command.c_str());
    
  }
  else if ( bam_count != 0 && fastq_count == 0 ){
    file_format = "bam";
    command 
      = DIR + "kmc -v -k" + std::to_string(kl) + " -t"+ std::to_string ( num_threads ) + " -m" + std::to_string ( maxRAM ) 
      + " -ci2 " + "-cx" + std::to_string(maxkfreq) + " -f" + file_format + " @" + infile_list + " sample " + tmp_dir + " 2>&1";
    std::cout << command << "\n" ; 
    print_file ( infile_list );
    system(command.c_str());
    
  }
  else if ( fastq_count != 0 && bam_count != 0){

    std::ofstream fout;
    fout.open("sample_fastq_file_list.txt");
    for ( auto i : fastq_files ){
      fout << i << "\n";
    }
    fout.close();
    fastq_files.clear();

    fout.open("sample_bam_file_list.txt");
    for ( auto i : bam_files ){
      fout << i << "\n";
    }
    fout.close();

    file_format = "q";
    command 
      = DIR + "kmc -v -k" + std::to_string(kl) + " -t"+ std::to_string ( num_threads ) + " -m" + std::to_string ( maxRAM ) 
      + " -ci2 -f" + file_format + " @sample_fastq_file_list.txt sample_fastq " + tmp_dir + " 2>&1";
    std::cout << command << "\n" ; 
    print_file ( "sample_fastq_file_list.txt" );
    system(command.c_str());
    
	
    file_format = "bam";
    command 
      = DIR + "kmc -v -k" + std::to_string(kl) + " -t"+ std::to_string ( num_threads ) + " -m" + std::to_string ( maxRAM ) 
      + " -ci2 -f" + file_format + " @sample_bam_file_list.txt sample_bam " + tmp_dir + " 2>&1";
    std::cout << command << "\n" ; 
    print_file ( "sample_bam_file_list.txt" );
    system(command.c_str());
    

    command = DIR + "kmc_tools simple sample_fastq sample_bam union sample 2>&1";
    std::cout << command << "\n" ; 
    system(command.c_str());
    if ( store_kmc == 0 ) {
      std::cout << "remove(\"sample_fastq.kmc_pre\")\n";
      remove("sample_fastq.kmc_pre");
      std::cout << "remove(\"sample_fastq.kmc_suf\")\n";
      remove("sample_fastq.kmc_suf");
      std::cout << "remove(\"sample_bam.kmc_pre\")\n";
      remove("sample_bam.kmc_pre");
      std::cout << "remove(\"sample_bam.kmc_suf\")\n";
      remove("sample_bam.kmc_suf");
    }
  }
  
  
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

    command = DIR + "kmc_tools transform sample histogram " + sample_hist + " 2>&1";
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
  std::cout << "k-mer cutoff: " << cutoff << "\n";


  ///////////////////////////////////////////////////////////////////////////////////////
  //
  // Applying filter
  //
  //std::cout << "Appying filter\n";
  command = DIR + "kmc_tools simple sample -ci" + std::to_string(cutoff) + " filter -ci1 kmers_subtract " + prefix + " 2>&1";
  std::cout << command << "\n" ;
  system(command.c_str() );
  if ( store_kmc == 0 ) {
    // std::cout << "remove(\"sample.kmc_pre\")\n";
    // remove("sample.kmc_pre");
    // std::cout << "remove(\"sample.kmc_suf\")\n";
    // remove("sample.kmc_suf");
    std::cout << "remove(\"filter.kmc_pre\")\n";
    remove("filter.kmc_pre");
    std::cout << "remove(\"filter.kmc_suf\")\n";
    remove("filter.kmc_suf");
  }

  ///////////////////////////////////////////////////////////////////////////////////////
  //
  // Print filter
  //
  //std::cout << "Writing sample-specific k-mers\n";
  command = DIR + "kmc_dump " + prefix + " " + prefix + ".txt" ;
  std::cout << command << "\n"; 
  system(command.c_str() );
  std::string kmc_pre=prefix+".kmc_pre";
  std::string kmc_suf=prefix+".kmc_suf";
  if ( store_kmc == 0 ) {
    std::cout << "remove(kmc_pre.c_str())\n";
    remove(kmc_pre.c_str());
    std::cout << "remove(kmc_suf.c_str())\n";
    remove(kmc_suf.c_str());
  }
  
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
