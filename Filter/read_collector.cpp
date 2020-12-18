//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#include "read_collector.hpp"


void read_collector_usage(){
  std::cout << "Usage: read_collector [options]" << "\n"
	    << "\n"
	    << "Required:\n"
	    << "\t" << "-f <string>\t" << "k-mer table\n"
	    << "\t" << "-1 <string>\t" << "1st fastq (gz supported)\n"
	    << "\t" << "-2 <string>\t" << "2nd fastq (gz supported)\n"
	    << "\t" << "           \t" << "-1 and -2 must be used simultaneously.\n"
	    << "\t" << "-b <string>\t" << "Bam file of paired-end reads\n"
	    << "\t" << "           \t" << "Do not use this option with -1 or -2.\n"
	    << "\n"
	    << "Optional:\n"
	    << "\t" << "-p <string>\t" << "Prefix of output files [filtered_read]\n"
	    << "\t" << "-t <int>   \t" << "Number of threads [8]\n"
	    << "\t" << "-l <int>   \t" << "K-mer size (<=32) [31]\n"
	    << "\n"
	    << "Contact:\n\tJang-il Sohn (sohnjangil@gmail.com)\n\tJin-Wu Nam (jwnam@hanyang.ac.kr)\n";
}

int main (int argc , char ** argv){

  if ( argc == 1 ){
    read_collector_usage();
    return 0;
  }

  int opt ;

  std::string input_1;
  std::string input_2;
  std::string input_b;
  int num_threads = 8;
  std::string prefix = "filtered_read";
  std::string kmer_table;
  int kl=31;

  while ( (opt = getopt ( argc, argv, "1:2:b:p:t:l:f:h" ) ) != -1 ){
    switch ( opt ) {
    case 'f': kmer_table = optarg; break;
    case '1': input_1 = optarg; break;
    case '2': input_2 = optarg; break;
    case 'b': input_b = optarg; break;
    case 'p': prefix = optarg; break;
    case 'l': kl = atoi(optarg); break;
    case 't': num_threads = atoi(optarg); break;
    case 'h': read_collector_usage(); return 0 ;
    default: std::cout << "\tInvalid option: " << optarg << "\n" ; read_collector_usage(); return 1;
    }
  }

  if ( kl > 32 || kl < 1){
    std::cout << "ERROR!!! -l must be positive number <=32\n";
    std::cout << "\n";
    read_collector_usage();
    return -1;
  }

  if ( kmer_table.size() == 0 ){
    std::cout << "ERROR!!! -f is required.\n";
    std::cout << "\n";
    read_collector_usage();
    return -1;
  }

  std::ifstream fin ( kmer_table.c_str() );
  if ( ! fin ){
    std::cout << "ERROR!!! There is no file: " << kmer_table << "\n";
    return -1 ;
  }
  fin.close();

  if ( input_1.size() != 0 && input_1.size() != 0 && input_b.size() != 0 ){
    std::cout << "WARNING!!! You input fastq and bam files simultaneously. Bam file " + input_b + " will be skipped.\n";
  }

  if ( input_1.size() != 0 && input_2.size() != 0 ){

    fin.open ( input_1.c_str() );
    if ( ! fin ){
      std::cout << "ERROR!!! There is no file: " << input_1.c_str() << "\n";
      return -1 ;
    }
    fin.close();
    
    fin.open ( input_2.c_str() );
    if ( ! fin ){
      std::cout << "ERROR!!! There is no file: " << input_2.c_str() << "\n";
      return -1 ;
    }
    fin.close();
    

    int gz_check = 0 ;
    if ( input_1.substr(input_1.size()-3) == ".gz" && input_2.substr(input_2.size()-3) == ".gz" ) {
      gz_check = 1;
    }

    if ( input_1.substr(input_1.size()-3) == ".gz" && input_2.substr(input_2.size()-3) != ".gz" ) {
      std::cout << "ERROR!!! Check input file names.\n\n";
      read_collector_usage();
      return -1;
    }

    if ( input_1.substr(input_1.size()-3) != ".gz" && input_2.substr(input_2.size()-3) == ".gz" ) {
      std::cout << "ERROR!!! Check input file names.\n\n";
      read_collector_usage();
      return -1;
    }
    collector (kmer_table, kl, prefix, input_1, input_2, gz_check, num_threads);
  }
  else if ( input_b.size() != 0 ){
    collector_bam (kmer_table, kl, prefix, input_b, num_threads);
  }
  else {
    std::cout << "ERROR!!! You missed input files.\n";
    std::cout << "\n";
    read_collector_usage();
    return -1;
  }
  
  return 0;
}
