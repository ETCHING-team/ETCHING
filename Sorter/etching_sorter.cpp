//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#include "etching_sorter.hpp"
#include <stdlib.h>
#include <random>

void usage(){
  std::cout
    << "Usage:\tetching_sorter [options] -i input.vcf -o output_prefix [options]\n"
    << "\n"
    << "Required:\n"
    << "\t" << "-i (string)\t" << "input vcf file (required)\n"
    << "\t" << "-o (string)\t" << "Prefix of output vcf file (required)\n"
    << "\n"
    << "Cut-off option:\n"
    << "\t" << "-c (double)\t" << "Cut-off parameter [0.4]\n"
    << "\n"
    << "Algorithm option:\n"
    << "\t" << "-R         \t" << "Random Forest [default]\n"
    << "\t" << "-X         \t" << "XGBoost\n"
    << "\t" << "-m         \t" << "Path to machine learning model [None]\n"
    << "\t" << "           \t" << "Use this if you want to specify machine learning model.\n"
    << "\n"
    << "I/O option:\n"
    << "\t" << "-T         \t" << "Tag SVs instead of removing low quality SVs [NONE]\n"
    << "\t" << "           \t" << "If you use this option, high confidential SVs will be tagged with\n"
    << "\t" << "           \t" << "\"PASS\", and low quality SVs will be tagged with \"LOWQUAL\" in\n"
    << "\t" << "           \t" << "FILTER field.\n"
    << "\n";
}
  
int main ( int argc , char ** argv ){
  if ( argc == 1 ){
    usage();
    return 0;
  }

  int opt;

  std::string infile;
  std::string prefix;
  std::string method;
  std::string path;
  int tagging=0;

  double alpha=-1;

  std::size_t sz;

  while ( (opt = getopt ( argc, argv, "i:o:c:m:RXT" ) ) != -1 ){
    switch ( opt ) {
    case 'i': infile=optarg; break; // infile name
    case 'o': prefix=optarg; break; // output prefix
    case 'c': alpha=std::stod(optarg,&sz); break; // cutoff
    case 'm': path=optarg; break; // path to machine learning model
    case 'R': method+="RandomForest"; break; // Random Forest
    case 'X': method+="XGBoost"; break; // XGBoost
    case 'T': tagging=1; break; // Tagging SVs instead of ramoving low quality SVs.
    default: std::cerr << "ERROR!!! Check options!!!\n\n" ; usage(); return 0 ; 
    }
  }
  
  if ( infile.size() == 0 ){
    std::cerr << "ERROR!!! -i (input file) is required\n"; 
    std::cerr << "------------------------------------\n";
    std::cerr << "\n";
    usage();
    return 0;
  }

  std::ifstream fin ( infile.c_str() );

  if ( ! fin ){
    std::cerr << "ERROR!!! There is no input file: " << infile << "\n";
    std::cerr << "-------------------------------------------------\n";
    std::cerr << "\n";
    usage();
    return 0;
  }

  fin.close();

  if ( method.size() == 0 ){
    method = "RandomForest";
    //method = "XGBoost";
  }
  else if (method == "RandomForestXGBoost" || method == "XGBoostRandomForest") {
    std::cerr << "ERROR!!! -R (Random Forest) and -X (XGBoost) can not be used at the same time\n";
    std::cerr << "-----------------------------------------------------------------------------\n";
    std::cerr << "\n";
    usage();
    return 0;
  }

  if ( alpha < 0 ){
    alpha = 0.4;
  }

  if ( prefix.size() != 0 ) prefix += ".";

  if ( path.size()== 0 ) {
    path = "${ETCHING_ML_PATH}";
  }

  std::string outfile = prefix + "etching_sorter.vcf";
  std::string feature_file = prefix + "feature_table.txt";
  std::string score_file = prefix + "score.txt";

  std::cout << "Method\t" << method << "\n";
  std::cout << "Cut-off\t" << alpha << "\n";
  std::cout << "Out-file\t" << outfile << "\n";
  std::string echo="echo -e \"Path-to-model\\t" + path + "\"\n";
  system ( echo.c_str() );
  //std::cout << "Path-to-model\t" << path << "\n";

  std::string command;
  std::random_device rd;
  std::uniform_int_distribution<int> dis(0, 999999);
  std::mt19937 gen(rd());
  std::string tmperr = "tmperr" + std::to_string ( dis(gen) );

  std::cout << "Calculating features\n";
  calc_feature(infile,feature_file);
  std::cout << "Scoring:\n";
  command = "scorer_" + method + " " + feature_file + " " + score_file + " " + path + " 2> " + tmperr + " ; rm " + tmperr;
  //std::cout << command << "\n";
  echo="echo "+command;
  system ( echo.c_str() );
  system ( command.c_str() );
  std::cout << "Removing false positives\n";
  razor ( infile, score_file, outfile, alpha, method, tagging);

  return 0;
}


