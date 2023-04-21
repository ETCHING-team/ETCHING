//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#include <cctype>
#include "read_collector.hpp"


void read_collector_usage(){
  std::cout << "Usage: read_collector [options]" << "\n"
	    << "\n"
	    << "Required:\n"
	    << "\t" << "-f <string>\t" << "k-mer counting table file consisting of \"k-mer\" and \"number\"\n"
	    << "\t" << "           \t" << "without a header\n"
	    << "\t" << "-1 <string>\t" << "1st fastq (gz supported)\n"
	    << "\n"
	    << "Optional:\n"
	    << "\t" << "-2 <string>\t" << "2nd fastq (gz supported)\n"
	    << "\t" << "-p <string>\t" << "Prefix of output files [filtered_read]\n"
	    << "\t" << "-t <int>   \t" << "Number of threads [16]\n"
	    << "\n"
	    << "Contact:\n\tJang-il Sohn (sohnjangil@gmail.com)\n\tJin-Wu Nam (jwnam@hanyang.ac.kr)\n";
}


int count_column(std::string infile){
  int  numberOfColumns=1;
  bool previousWasSpace=false;

  std::string line;

  std::ifstream fin ( infile.c_str() );
  std::getline ( fin , line );
  fin.close();


  for(std::size_t i=0; i<line.size(); i++){
    if(line[i] == ' ' || line[i] == '\t'){
      if(!previousWasSpace)
	numberOfColumns++;

      previousWasSpace = true;
    } else {
      previousWasSpace = false;
    }   
  }   

  return numberOfColumns;
}


int get_kmer_size ( std::string kmer_table_fname ) {
  std::string kmer;
  std::ifstream fin ( kmer_table_fname.c_str() );
  if ( ! fin ){
    std::cerr << "ERROR!!! read_collector: There is no file: " << kmer_table_fname << "\n";
    return 1 ;
  }
  fin >> kmer;
  int kl = kmer.size();
  if ( kmer.find_first_not_of("ACTGacgtURYSWKMBDHVN.") != std::string::npos ){
    std::cerr << "ERROR!!! read_collector: Improper k-mer table: " << kmer_table_fname << std::endl;
    exit(1);
  }
  if ( kl > 32 ) {
    std::cerr << "ERROR!!! read_collector: k-mer size must be <=32. Please check your k-mer table file: " << kmer_table_fname << std::endl;
    exit(1);
  }
  fin.close();


  if ( count_column ( kmer_table_fname ) != 2 ){
    std::cerr << "ERROR!!! read_collector: Please check k-mer table file: " << kmer_table_fname << "\n";
    std::cerr << "It must be two-column of k-mer and count: " << std::endl;
    std::cerr << "Ex)" << std::endl;
    std::cerr << "AATGGTGAACTAGAAGTCCA 12" << std::endl;
    std::cerr << "AATGGTGAACTAGAAGTCCC 5" << std::endl;
    std::cerr << "AATGGTGAACTAGAAGTCGT 11" << std::endl;
    std::cerr << "." << std::endl;
    std::cerr << "." << std::endl;
    std::cerr << "." << std::endl;
    exit(1);
  }


  return kl;
}


int main (int argc , char ** argv){

  if ( argc == 1 ){
    read_collector_usage();
    return 1;
  }

  int opt ;

  std::string input_1;
  std::string input_2;
  int num_threads = 16;
  std::string prefix = "filtered_read";
  std::string kmer_table;
  std::string dummy;
  while ( (opt = getopt ( argc, argv, "1:2:p:t:f:l:h" ) ) != -1 ){
    switch ( opt ) {
    case 'f': kmer_table = optarg; break;
    case '1': input_1 = optarg; break;
    case '2': input_2 = optarg; break;
    case 'p': prefix = optarg; break;
    case 't': num_threads = atoi(optarg); break;
    case 'l': dummy = optarg; break;
    case 'h': read_collector_usage(); return 1 ;
    default: std::cerr << "\tInvalid option: " << optarg << "\n" ; read_collector_usage(); return 1;
    }
  }

  if ( kmer_table.size() == 0 ){
    std::cerr << "ERROR!!! read_collector: -f is required.\n";
    std::cerr << "\n";
    read_collector_usage();
    return 1;
  }


  int kl = get_kmer_size(kmer_table);
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Pair mode
  //

  if ( input_1.size() != 0 && input_2.size() != 0 ){
    std::ifstream fin;
    fin.open ( input_1.c_str() );
    if ( ! fin ){
      std::cerr << "ERROR!!! read_collector: There is no file: " << input_1.c_str() << "\n";
      return 1 ;
    }
    fin.close();
    
    fin.open ( input_2.c_str() );
    if ( ! fin ){
      std::cerr << "ERROR!!! read_collector: There is no file: " << input_2.c_str() << "\n";
      return 1 ;
    }
    fin.close();
    

    int gz_check = 0 ;
    if ( input_1.substr(input_1.size()-3) == ".gz" && input_2.substr(input_2.size()-3) == ".gz" ) {
      gz_check = 1;
    }

    if ( input_1.substr(input_1.size()-3) == ".gz" && input_2.substr(input_2.size()-3) != ".gz" ) {
      std::cerr << "ERROR!!! read_collector: Check input file names.\n\n";
      read_collector_usage();
      return 1;
    }

    if ( input_1.substr(input_1.size()-3) != ".gz" && input_2.substr(input_2.size()-3) == ".gz" ) {
      std::cerr << "ERROR!!! read_collector: Check input file names.\n\n";
      read_collector_usage();
      return 1;
    }

    // Run collector
    int check_collector = collector (kmer_table, kl, prefix, input_1, input_2, gz_check, num_threads);
    if ( check_collector != 0 ){
      std::cerr << "ETCHING-Filter was abnormally finished at the collector function in read_collector\n";
      return 1;
    }
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // single mode
  //
  else if ( input_1.size() != 0 && input_2.size() == 0 ){
    std::ifstream fin;
    fin.open ( input_1.c_str() );
    if ( ! fin ){
      std::cerr << "ERROR!!! read_collector: There is no file: " << input_1.c_str() << "\n";
      return 1;
    }
    fin.close();
    
    int gz_check = 0 ;
    if ( input_1.substr(input_1.size()-3) == ".gz" ){
      gz_check = 1;
    }

    // Run collector
    // k-mer table mode
    int check_collector_single = collector_single (kmer_table, kl, prefix, input_1, gz_check, num_threads);
    if ( check_collector_single != 0 ){
      std::cerr << "ETCHING-Filter was abnormally finished at the collector function in read_collector\n";
      return 1;
    }
  }
  
  return 0;
}
