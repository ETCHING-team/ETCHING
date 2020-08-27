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
	    << "\n"
	    << "Optional:\n"
	    << "\t" << "-o <string>\t" << "Prefix of output files [filtered_read]\n"
	    << "\t" << "-t <int>   \t" << "Number of threads [8]\n"
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
  int num_threads = 8;
  std::string prefix = "filtered_read";
  std::string kmer_table;
  int kl=31;

  while ( (opt = getopt ( argc, argv, "1:2:o:t:f:h" ) ) != -1 ){
    switch ( opt ) {
    case 'f': kmer_table = optarg; break;
    case '1': input_1 = optarg; break;
    case '2': input_2 = optarg; break;
    case 'o': prefix = optarg; break;
    case 't': num_threads = atoi(optarg); break;
    case 'h': read_collector_usage(); return 0 ;
    default: std::cerr << "\tInvalid option: " << optarg << "\n" ; read_collector_usage(); return 1;
    }
  }

  if ( kmer_table.size() == 0 ){
    std::cerr << "ERROR!!! -f is required.\n";
    std::cerr << "\n";
    read_collector_usage();
    return 0;
  }

  if (input_1.size()==0 || input_2.size()==0){
    std::cerr << "ERROR!!! -1 and -2 are required.\n";
    std::cerr << "\n";
    read_collector_usage();
    return 0;
  }

  int gz_check = 0 ;
  if ( input_1.substr(input_1.size()-3) == ".gz" && input_2.substr(input_2.size()-3) == ".gz" ) gz_check = 1;
  if ( input_1.substr(input_1.size()-3) == ".gz" && input_2.substr(input_2.size()-3) != ".gz" ) {
    std::cerr << "ERROR!!! Check input file names.\n\n";
    read_collector_usage();
    return 0;
  }
  if ( input_1.substr(input_1.size()-3) != ".gz" && input_2.substr(input_2.size()-3) == ".gz" ) {
    std::cerr << "ERROR!!! Check input file names.\n\n";
    read_collector_usage();
    return 0;
  }

  collector (kmer_table, kl, prefix, input_1, input_2, gz_check, num_threads);

}
