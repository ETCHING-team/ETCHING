//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------


#include "etching_caller.hpp"
#include <chrono>
#include <ctime>

#include "CPU_TIME_measure.hpp"
#include "Peak_memory_measure.hpp"

int main(int argc, char ** argv){

  if(argc == 1 ){
    caller_usage();
    return 0 ;
  }

  double cputime0 = get_cpu_time();
  double walltime0 = get_wall_time();

  int opt ;

  std::string infile_pref; // for step 2 and 3
  std::string input_bam;
  // std::string normal_bam;
  std::string prefix;
  std::string genome;
 
  // int num_threads = 8;
  
  double purity = 0.75;
  double seqdep = 0 ;

  int insert_size = 500; // -I insert_size
  bool typing = 1 ;
  bool scanall = 1 ; 
  bool rescue = 0 ;

  int step=1;

  std::string read_orientation = "FR";
  // std::string data_type; // WGS, WTS, or PANEL
  
  std::size_t sz;

  std::string confile;


  // while ( (opt = getopt ( argc, argv, "b:o:g:t:BS:i:D:P:T:I:O:c:r:a:f:x:q:e:C:R:A:F:X:Q:E:h" ) ) != -1 ){
  while ( (opt = getopt ( argc, argv, "b:o:g:BS:i:ARD:P:I:O:c:h" ) ) != -1 ){
    switch ( opt ) {
    case 'b': input_bam=optarg; break;
    case 'o': prefix=optarg; break;
    case 'g': genome=optarg; break;
    case 'B': typing = 0; break;
    case 'S': step = atoi(optarg); if ( step < 1 || step > 3) { caller_usage(); return 0; } break;
    case 'i': infile_pref = optarg; break;
    case 'A': scanall = 1; break;
    case 'R': rescue = 1; break;
    case 'D': seqdep = std::stod(optarg,&sz); break;
    case 'P': purity = std::stod(optarg,&sz); break;
    case 'I': insert_size = atoi(optarg); break;
    case 'O': read_orientation = optarg; break;
    case 'c': confile = optarg; break;
    case 'h': caller_usage(); return 0 ;

    default: std::cout << "\tInvalid option\n"; caller_usage(); return 1;

    }
  }

  if ( input_bam.size() == 0 ){
    std::cout << "------------------------------------------------------------------------------\n";
    std::cout << "ERROR!!! In option -i : Input bam file is required.\n";
    std::cout << "------------------------------------------------------------------------------\n\n";
    caller_usage();
    return 0;
  }
  if ( prefix.size() == 0 ){
    std::cout << "------------------------------------------------------------------------------\n";
    std::cout << "ERROR!!! In option -o : Prefix of output file is required.\n";
    std::cout << "------------------------------------------------------------------------------\n\n";
    caller_usage();
    return 0;
  }
  if ( genome.size() == 0 ){
    std::cout << "------------------------------------------------------------------------------\n";
    std::cout << "ERROR!!! In option -g : Reference genome file is required.\n";
    std::cout << "------------------------------------------------------------------------------\n\n";
    caller_usage();
    return 0;
  }
  if ( read_orientation != "FR" && read_orientation != "RF" ) {
    std::cout << "------------------------------------------------------------------------------\n";
    std::cout << "ERROR!!! In option -O : Read-orientation must be one of FR or RF.\n";
    std::cout << "------------------------------------------------------------------------------\n\n";
    caller_usage();
    return 0;
  }
  if ( purity <= 0 || purity > 1 ){
    std::cout << "------------------------------------------------------------------------------\n";
    std::cout << "ERROR!!! In option -P : Tumor purity must be from 0 to 1.\n";
    std::cout << "------------------------------------------------------------------------------\n\n";
    caller_usage();
    return 0;
  }
  if ( seqdep < 0 ){
    std::cout << "------------------------------------------------------------------------------\n";
    std::cout << "ERROR!!! In option -C : Sequencing coverage (or depth) must be greater then 0.\n";
    std::cout << "------------------------------------------------------------------------------\n\n";
    caller_usage();
    return 0;
  }
  
  std::ifstream fin ( input_bam.c_str() );
  if ( !fin.good() ) {
    std::cout << "------------------------------------------------------------------------------\n";
    std::cout << "ERROR!!! In option -b : No input bam file : " << input_bam << "\n";
    std::cout << "------------------------------------------------------------------------------\n\n";
    caller_usage();
    return 0;
  }
  fin.close();

  fin.open( genome.c_str() );
  if ( !fin.good() ) {
    std::cout << "------------------------------------------------------------------------------\n";
    std::cout << "ERROR!!! In option -g : No reference genome file : " << genome << "\n";
    std::cout << "------------------------------------------------------------------------------\n\n";
    caller_usage();
    return 0;
  }
  fin.close();

  //////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Program start
  //

  int Step = step;


  //////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Printing parameters
  //
  auto start = std::chrono::system_clock::now();
  std::time_t start_time = std::chrono::system_clock::to_time_t(start);
  
  std::cout << "\n";
  std::cout << "===============================================================\n";
  std::cout << "[Start]\n" << std::ctime(&start_time) << "\n";
  std::cout << "[Required parameters]\n";
  std::cout << "Input bam file:       -b " << input_bam << "\n";
  std::cout << "Prefix of output:     -o " << prefix << "\n";
  std::cout << "Reference genome:     -g " << genome << "\n";
  std::cout << "Insert-size:          -I " << insert_size << "\n";
  std::cout << "Read-orientation:     -O " << read_orientation << "\n";
  if ( typing == 0 ) std::cout << "BND only:             -B " << "Print only BND" << "\n";
  //if ( scanall ) std::cout << "Scan all split-read:  -A " << "\n";
  if ( rescue )  std::cout << "Rescue SVs:           -R " << "\n";
  std::cout << "\n";
  std::cout << "===============================================================\n";

  //////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Step 1
  //

  if ( Step ==1 ){
    std::string bp_file;
    find_path ( input_bam, prefix, insert_size, scanall, rescue);
    Step ++ ;

  }

  std::cout << "===============================================================\n";
  //////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Step 2
  //

  if ( Step == 2){
    if (step == 1 ) {
      infile_pref = prefix ;
    }
    bp_to_vcf ( genome, input_bam, infile_pref,	prefix,	insert_size, typing, read_orientation, purity, seqdep);

    Step ++;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////
  //
  // End of program
  //

  double walltime1 = get_wall_time();
  double cputime1 = get_cpu_time();

  auto end = std::chrono::system_clock::now();
  
  std::cout << "===============================================================\n";
  std::cout << "Total CPU time: " << cputime1 - cputime0 << " sec\n";
  std::cout << "Total wall-clock time: " << walltime1 - walltime0 << " sec\n";
  std::cout << "Peak memory: " << std::to_string((double)get_peak_memory()/(1024*1024)) << " GB\n";
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  std::cout << "[End]\n" << std::ctime(&end_time) << "\n";
  std::cout << "[Output files]\n";
  std::cout << prefix << ".BND.vcf\n";
  if ( typing ) std::cout << prefix << ".SV.vcf\n";
  std::cout << "===============================================================\n";

  return 0;
}
