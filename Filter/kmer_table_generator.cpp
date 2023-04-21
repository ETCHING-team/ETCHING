//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#ifndef ETCHING_KMER_TABLE_GENERATOR
#define ETCHING_KMER_TABLE_GENERATOR

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include <unistd.h>
#include <sys/stat.h>

#include "CPU_TIME_measure.hpp"

void filter_usage(){
  std::cout << "Usage: kmer_table_generator [options]\n"
	    << "\n"
	    << "Options:\n"
	    << "\t" << "-o <string>\t" << "Prefix of final output file name [required]\n"
	    << "\t" << "           \t" << "\tPREFIX.kmer_table.txt\n"
	    << "\t" << "-i <string>\t" << "List of sample (or tumor) in fastq(.gz)\n"
	    << "\t" << "           \t" << "Uncomfortable with -I\n"
	    << "\t" << "-I <string>\t" << "(For reuse) Prefix of existing sample k-mer set:\n"
	    << "\t" << "           \t" << "\tPREFIX.sample.kmc_pre/suf\n"
	    << "\t" << "           \t" << "Uncomfortable with -i\n"
	    << "\t" << "-c <string>\t" << "List of control (or matched normal) in fastq(.gz)\n"
	    << "\t" << "           \t" << "Uncomfortable with -C\n"
	    << "\t" << "-C <string>\t" << "(For reuse) Prefix of existing control k-mer set:\n"
	    << "\t" << "           \t" << "\tPREFIX.control.kmc_pre/suf\n"
	    << "\t" << "           \t" << "Uncomfortable with -c\n"
	    << "\t" << "-a <string>\t" << "Add k-mer set. To add PGK2, you must have\n"
	    << "\t" << "           \t" << "PGK2.kmc_pre/suf\n"
	    << "\t" << "-l <string>\t" << "k-mer size [31]\n"
	    << "\t" << "-K <int>   \t" << "k-mer counting cut-off for remove sequencing error [3]\n"
	    << "\t" << "-S <int>   \t" << "Maximum counting value [255]\n"
	    << "\t" << "-M <int>   \t" << "The k-mer exceeding this will be ignored [10000]\n"
	    << "\t" << "-m <int>   \t" << "Maximum RAM memory in GB [12]\n"
	    << "\t" << "-t <int>   \t" << "Number of thread [8]\n"
	    << "\n"
	    << "Contact:\n\tJang-il Sohn (sohnjangil@gmail.com)\n\tJin-Wu Nam (jwnam@hanyang.ac.kr)\n";
}

inline bool check_exist(const std::string& name) {
  struct stat buffer;   
  return(stat(name.c_str(), &buffer) == 0); 
}



int main( int argc , char ** argv ){
  if( argc == 1 ){
    filter_usage();
    return 0;
  }

  std::string sample_list_file;
  std::string control_list_file;
  std::string prefix="filtered_kmer";

  std::string given_sample;
  std::string given_control;
  std::string additional;

  std::string result;

  int64_t maxRAM = 12;
  int64_t num_threads = 8;

  int64_t cutoff = 3;
  int64_t maxcount = 255;
  int64_t maxkfreq = 10000;
  int64_t kl = 31;

  int64_t opt = 0;
  while((opt = getopt( argc , argv , "i:I:p:o:c:C:a:m:t:K:S:M:l:D:h") ) != -1 ){
    switch(opt){
    case 'i': sample_list_file = optarg; break;
    case 'I': given_sample = optarg; break;
    case 'p': prefix = optarg; break;
    case 'o': prefix = optarg; break;
    case 'c': control_list_file = optarg; break;
    case 'C': given_control = optarg; break;
    case 'a': additional = optarg; break;
    case 'm': maxRAM = atoi( optarg ); break;
    case 't': num_threads = atoi( optarg ); break;
    case 'K': cutoff = atoi( optarg ); break;
    case 'S': maxcount = atoi( optarg ); break;
    case 'M': maxkfreq = atoi( optarg ); break;
    case 'l': kl = atoi(optarg ); break;
    case 'h': filter_usage(); return 0;
    default: std::cerr << "ERROR!!! Invalid option: " << optarg << "\n"; return 1;
    }
  }

  if( system("kmc > /dev/null 2>&1") != 0 ){
    std::cerr << "ERROR!!! We cannot find kmc in your PATH: ";
    system("echo $PATH");
    return 1;
  }
  if( system("kmc_tools > /dev/null 2>&1") != 0 ){
    std::cerr << "ERROR!!! We cannot find kmc_tools in your PATH: ";
    system("echo $PATH");
    return 1;
  }
  if( system("kmc_dump > /dev/null 2>&1") != 0 ){
    std::cerr << "ERROR!!! We cannot find kmc_dump in your PATH: ";
    system("echo $PATH");
    return 1;
  }

  ///////////////////////////////////////////////////////////////////////////////////////
  //
  // Option check
  //

  if( cutoff > maxkfreq ){
    std::cerr << "ERROR!!! -K must be less than -M\n";
    return 1;
  }

  if( cutoff > maxcount ){
    std::cerr << "ERROR!!! -K must be less than -S\n";
    return 1;
  }

  if( sample_list_file.size() + given_sample.size() == 0){
    std::cerr << "ERROR!!! -i or -I is required.\n";
    return 1;
  }

  if( control_list_file.size() + given_control.size() + additional.size() == 0){
    std::cerr << "ERROR!!! -c, -C, or -a were required\n";
    return 1;
  }

  std::string tmp;
  std::ifstream fin( sample_list_file );
  
  if( sample_list_file.size() != 0 ){
    std::cout << "## Check sample files:\n";
    while( fin >> tmp ){
      std::cout << "##\t" << tmp << std::endl;
      if( tmp.substr(tmp.size()-6) != ".fastq" && 
	  tmp.substr(tmp.size()-3) != ".fq" && 
	  tmp.substr(tmp.size()-9) != ".fastq.gz" && 
	  tmp.substr(tmp.size()-6) != ".fq.gz"){
	std::cerr << "ERROR!!! Please check files in " << tmp <<".\n";
	std::cerr << "This file must be fastq(.gz) format.\n";
      }
    }
    fin.close();
  }

  if( given_sample.size() != 0 ){
    tmp = given_sample + ".kmc_pre";
    std::ifstream fin;
    std::cout << "## Check sample k-mer set:" << std::endl;
    std::cout << "##\t" << tmp << std::endl;
    fin.open( tmp.c_str() );
    if( ! fin ){
      std::cerr << "ERROR!!! There is no " + given_sample + ".kmc_pre.\n";
      return 1;
    }
    fin.close();
    tmp = given_sample + ".kmc_suf";
    std::cout << "##\t" << tmp << std::endl;
    fin.open( tmp.c_str() );
    if( ! fin ){
      std::cerr << "ERROR!!! There is no " + given_sample + ".kmc_suf.\n";
      return 1;
    }
  }

  
  if( control_list_file.size() != 0 ){
    std::ifstream fin( control_list_file.c_str() );
    std::cout << "## Check control files:" << std::endl;
    while( fin >> tmp){
      std::cout << "##\t" << tmp << std::endl;
      if( tmp.substr(tmp.size()-6) != ".fastq" && 
	  tmp.substr(tmp.size()-3) != ".fq" && 
	  tmp.substr(tmp.size()-9) != ".fastq.gz" && 
	  tmp.substr(tmp.size()-6) != ".fq.gz"){
	std::cerr << "ERROR!!! Please check files in " << tmp <<".\n";
	std::cerr << "This file must be fastq(.gz) format.\n";
      }
    }
    fin.close();
  }
  
  if( given_control.size() != 0 ){
    tmp = given_control + ".kmc_pre";
    std::ifstream fin;
    std::cout << "## Check control k-mer set:" << std::endl;
    std::cout << "##\t" << tmp << std::endl;
    fin.open( tmp.c_str() );
    if( ! fin ){
      std::cerr << "ERROR!!! There is no " + given_control + ".kmc_pre.\n";
      return 1;
    }
    fin.close();
    tmp = given_control + ".kmc_suf";
    std::cout << "##\t" << tmp << std::endl;
    fin.open( tmp.c_str() );
    if( ! fin ){
      std::cerr << "ERROR!!! There is no " + given_control + ".kmc_suf.\n";
      return 1;
    }
  }


  if( additional.size() != 0 ){
    tmp = additional + ".kmc_pre";
    std::ifstream fin;
    std::cout << "## Check additional k-mer set:" << std::endl;
    std::cout << "##\t" << tmp << std::endl;
    fin.open( tmp.c_str() );
    if( ! fin ){
      std::cerr << "ERROR!!! There is no " + additional + ".kmc_pre.\n";
      return 1;
    }
    fin.close();
    tmp = additional + ".kmc_suf";
    std::cout << "##\t" << tmp << std::endl;
    fin.open( tmp.c_str() );
    if( ! fin ){
      std::cerr << "ERROR!!! There is no " + additional + ".kmc_suf.\n";
      return 1;
    }
  }

  std::cout << "## All input files are OK.\n\n";

  ///////////////////////////////////////////////////////////////////////////////////////
  //
  // Generate random tmp, filter, control, and sample names
  //

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::mt19937 generator(seed);

  std::string RN;

  std::string tmp_dir;

  std::string sample_prefix = prefix + ".sample";  // meaning tumor
  std::string control_prefix= prefix + ".control"; // meaning normal
  std::string result_prefix = prefix;

  // check if tmp_dir_RN
  // if one of them exists, try another RN
  int max_try_limit = 100;
  while( 1 ){
    RN = std::to_string(generator());

    tmp_dir = "tmp_" + RN;
    int bad_RN = check_exist( tmp_dir );

    if( bad_RN == 0 ) break;
    if( max_try_limit ) {
      std::cerr << "ERROR!!! " << __FILE__ << " (line " << __LINE__ << ")::exceed max_try_limit " << max_try_limit << "::Highly likely a bug." << std::endl;
      return 1;
    }
  }

  std::string command = "mkdir -p " + tmp_dir;
  std::cout << command << "\n"; 
  if ( system(command.c_str()) != 0 ){
    std::cerr << "Failed!!! mkdir" << std::endl;
    return 1;
  }

  ///////////////////////////////////////////////////////////////////////////////////////
  //
  // Basic KMC command
  //

  std::string KMC
    = "kmc -v -k" + std::to_string(kl)
    + " -t"+ std::to_string(num_threads)
    + " -m" + std::to_string(maxRAM)
    + " -cx" + std::to_string(maxkfreq)
    + " -ci2" + 
    + " -fq";
  ///////////////////////////////////////////////////////////////////////////////////////
  //
  // Making tumor k-mer table
  //
  
  if( sample_list_file.size() !=0 ){
    command 
      = KMC + " @" + sample_list_file + " " + sample_prefix + " " + tmp_dir + " 2>&1 | sed \"s/\\r/\\n/g\" > " + sample_prefix  + ".kmc.log";
    std::cout << command << "\n"; 
    if( system(command.c_str()) != 0 ){
      std::cerr << "Failed!!! KMC did not work properly." << std::endl;
      std::cerr << "You may check: " << sample_prefix << ".kmc.log\n";
      return 1;
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////
  //
  // Making normal filter
  //

  if( control_list_file.size() !=0 ){
    command
      = KMC + " @" + control_list_file + " " + control_prefix + " " + tmp_dir + " 2>&1 | sed \"s/\\r/\\n/g\" > " + control_prefix  + ".kmc.log";
    std::cout << command << "\n"; 
    if( system(command.c_str()) != 0 ){
      std::cerr << "Failed!!! kmc" << std::endl;
      std::cerr << "You may check: " << control_prefix << ".kmc.log\n";
      return 1;
    }
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  //
  // Make config file for kmc_tools
  //


  std::string config_fname = prefix + ".kmc.conf";

  std::string sample_config;
  std::string control_config;
  std::string additional_config;
  std::string result_config;

  if ( sample_list_file.size() != 0 && given_sample.size() == 0 ){
    sample_config 
      = "sample = " + sample_prefix
      + " -cx" + std::to_string(maxkfreq)
      + " -ci" + std::to_string(cutoff);
  }
  else if ( sample_list_file.size() == 0 && given_sample.size() != 0 ){
    sample_config 
      = "sample = " + given_sample
      + " -cx" + std::to_string(maxkfreq)
      + " -ci" + std::to_string(cutoff);
  }
  else {
    std::cerr << "ERROR!!! " << __FILE__ << "(line " << __LINE__ << "):: no sample given\n";
    return 1;
  }

  if ( control_list_file.size() != 0 && given_control.size() == 0 ){
    control_config 
      ="control = " + control_prefix
      + " -cx" + std::to_string(maxkfreq)
      + " -ci1";
  }
  else if ( control_list_file.size() == 0 && given_control.size() != 0 ){
    control_config 
      = "control = " + given_control
      + " -cx" + std::to_string(maxkfreq)
      + " -ci1";
  }

  if ( additional.size() != 0 ){
    additional_config = "additional = " + additional + " -ci1";
  }

  if ( control_config.size() != 0 && additional_config.size() != 0 ){
    result_config = result_prefix + " = sample - ( control + additional )";
  }
  else if ( control_config.size() == 0 && additional_config.size() != 0 ){
    result_config = result_prefix + " = sample - additional";
  }
  else if ( control_config.size() != 0 && additional_config.size() == 0 ){
    result_config = result_prefix + " = sample - control";
  }
  else{
    std::cerr << "ERROR!!! " << __FILE__ << " (line " << __LINE__ << "):: no control nor additional\n";
    return 1;
  }


  std::ofstream fout ( config_fname.c_str() );
  fout << "INPUT:\n";
  fout << sample_config << "\n";
  if ( control_config.size() != 0 ) fout << control_config << "\n";
  if ( additional_config.size() != 0 ) fout << additional_config << "\n";
  fout << "OUTPUT:\n";
  fout << result_config << "\n";

  fout.close();

  ///////////////////////////////////////////////////////////////////////////////////////
  //
  // Applying filter (kmc_tools complex)
  //

  command = "kmc_tools complex " + config_fname + " 2>&1 | sed \"s/\\r/\\n/g\" > " + prefix + ".kmc_tools.log";
  std::cout << command << "\n";
  if ( system(command.c_str()) != 0 ){
    std::cerr << "Failed!!! kmc_tools complex" << std::endl;
    std::cerr << "You may check: " << prefix << ".kmc_tools.log\n";
    return 1;    
  }
  ///////////////////////////////////////////////////////////////////////////////////////
  //
  // Print filter
  //
  command = "kmc_dump " + result_prefix + " " + prefix + ".kmer_table.txt ";
  std::cout << command << "\n"; 
  if ( system(command.c_str()) != 0 ){
    std::cerr << "Failed!!! kmc_dump" << std::endl;
    return 1;    
  }
  std::string kmc_pre=prefix+".kmc_pre";
  std::string kmc_suf=prefix+".kmc_suf";
  
  ///////////////////////////////////////////////////////////////////////////////////////
  //
  // End
  //

  command = "rm -rf " + tmp_dir;

  std::cout << command << "\n";
  if( system(command.c_str()) != 0 ){
    std::cerr << "Failed!!! remove tmp directory" << command << std::endl;
    return 1;
  }

  return 0;
}


#endif
