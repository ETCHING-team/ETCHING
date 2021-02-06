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
	    << "\t" << "-f <string>\t" << "k-mer table (or a list of read IDs to be selected)\n"
	    << "\t" << "-1 <string>\t" << "1st fastq (gz supported)\n"
	    << "\t" << "-2 <string>\t" << "2nd fastq (gz supported)\n"
	    << "\t" << "-b <string>\t" << "Bam file of paired-end reads\n"
	    << "\t" << "           \t" << "Do not use this option with -1 or -2.\n"
	    << "\n"
	    << "Optional:\n"
	    << "\t" << "-p <string>\t" << "Prefix of output files [filtered_read]\n"
	    << "\t" << "-t <int>   \t" << "Number of threads [8]\n"
	    << "\t" << "-l <int>   \t" << "K-mer size (<=32) [31]\n"
	    << "\t" << "           \t" << "-1 and -2 must be used simultaneously.\n"
	    << "\t" << "-F         \t" << "Fast-bam mode for -b option\n"
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
  bool fast_bam = 0 ;
  int kl=31;

  while ( (opt = getopt ( argc, argv, "1:2:b:p:t:l:f:Fh" ) ) != -1 ){
    switch ( opt ) {
    case 'f': kmer_table = optarg; break;
    case '1': input_1 = optarg; break;
    case '2': input_2 = optarg; break;
    case 'b': input_b = optarg; break;
    case 'p': prefix = optarg; break;
    case 'l': kl = atoi(optarg); break;
    case 't': num_threads = atoi(optarg); break;
    case 'F': fast_bam = 1; break;
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

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Pair mode
  //

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


    // Check if k-mer table mode or id mode
    bool id_mode = 0 ;
    if ( count_column ( kmer_table ) == 2 ){
      id_mode = 0;
    }
    else if ( count_column ( kmer_table ) == 1 ){
      id_mode = 1;
    }
    else {
      std::cout << "ERROR!!! Please check a file: " << kmer_table << "\n";
      return count_column ( kmer_table );
    }

    // Run collector
    // k-mer table mode
    if ( ! id_mode ){
      std::cout << "Entering k-mer table mode\n";
      int check_collector = collector (kmer_table, kl, prefix, input_1, input_2, gz_check, num_threads);
      if ( check_collector != 0 ){
	std::cout << "ETCHING-Filter was abnormally finished at the collector function in read_collector\n";
	return check_collector;
      }
    }
    // id table mode
    else{
      std::cout << "Entering ID table mode\n";
      std::string id_table = kmer_table;
      int check_collector_id_mode = collector_id_mode (id_table, prefix, input_1, input_2, gz_check, num_threads);
      if ( check_collector_id_mode != 0 ){
	std::cout << "ETCHING-Filter was abnormally finished at the collector function in read_collector\n";
	return check_collector_id_mode;
      }
    }
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // single mode
  //
  else if ( input_1.size() != 0 && input_2.size() == 0 ){ 
    std::cout << "Entering single mode\n";
    fin.open ( input_1.c_str() );
    if ( ! fin ){
      std::cout << "ERROR!!! There is no file: " << input_1.c_str() << "\n";
      return -1 ;
    }
    fin.close();
    
    int gz_check = 0 ;
    if ( input_1.substr(input_1.size()-3) == ".gz" ){
      gz_check = 1;
    }

    // Check if k-mer table mode or id mode
    bool id_mode = 0 ;

    if ( count_column ( kmer_table ) == 2 ){
      id_mode = 0;
    }
    else if ( count_column ( kmer_table ) == 1 ){
      id_mode = 1;
    }
    else {
      std::cout << "ERROR!!! Please check a file: " << kmer_table << "\n";
      return count_column ( kmer_table );
    }

    // Run collector
    // k-mer table mode
    if ( ! id_mode ){
      int check_collector_single = collector_single (kmer_table, kl, prefix, input_1, gz_check, num_threads);
      if ( check_collector_single != 0 ){
	std::cout << "ETCHING-Filter was abnormally finished at the collector function in read_collector\n";
	return check_collector_single;
      }
    }
    else{
      std::string id_table = kmer_table;
      int check_collector_id_mode_single = collector_id_mode_single (id_table, prefix, input_1, gz_check, num_threads);
      if ( check_collector_id_mode_single != 0 ){
	std::cout << "ETCHING-Filter was abnormally finished at the collector function in read_collector\n";
	return check_collector_id_mode_single;
      }
    }
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // bam mode
  //
  else if ( input_b.size() != 0 ){
    int check_collector = collector_bam (kmer_table, kl, prefix, input_b, num_threads, fast_bam);
    if ( check_collector != 0 ){
      std::cout << "ETCHING-Filter was abnormally finished at the collector_bam function in read_collector\n";
      return check_collector;
    }
  }
  else {
    std::cout << "ERROR!!! You missed input files.\n";
    std::cout << "\n";
    read_collector_usage();
    return -1;
  }
  
  return 0;
}
