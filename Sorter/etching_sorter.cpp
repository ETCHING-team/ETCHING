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
    << "Options:\n"
    << "\t" << "-c (double)\t" << "Cut-off parameter [0.4]\n"
    << "\t" << "-R         \t" << "Random Forest [default]\n"
    << "\t" << "-X         \t" << "XGBoost\n"
    << "\t" << "-m         \t" << "Path to machine learning mode\n"
    // << "\n"
    // << "I/O option:\n"
    // << "\t" << "-Q         \t" << "Tag SVs with \"PASS\" or \"LOWQUAL\" instead of removing low quality SVs [NONE]\n"
    << "\n";
}


int python_version(std::string prefix){
  std::string fname = prefix + "python_version";
  std::string command = "python3 --version > " + fname + " 2>&1";
  system ( command.c_str() );
  std::ifstream fin ( fname.c_str() );
  std::string tmp;
  fin >> tmp >> tmp;
  tmp = tmp[0];
  return atoi ( tmp.c_str() );
}


std::string python_module_check(std::string prefix, std::string method){
  std::string output;

  std::vector < std::string > pm;
  std::vector < std::string > module;

  pm.push_back("import pickle");
  pm.push_back("import pandas");
  pm.push_back("import numpy");
  pm.push_back("from sklearn import metrics");
  if ( method == "RandomForest")
    pm.push_back("import skranger");
  if ( method == "XGBoost")
    pm.push_back("import xgboost");

  module.push_back("pickle");
  module.push_back("pandas");
  module.push_back("numpy");
  module.push_back("sklearn");
  if ( method == "RandomForest")
    module.push_back("skranger");
  if ( method == "XGBoost")
    module.push_back("xgboost");

  for ( std::size_t i = 0 ; i < pm.size() ; i ++ ){
    std::string fname = prefix + "python_module_check";
    std::string command = "python -c '" + pm[i] + "' 2>&1 | wc -l > " + fname;
    system ( command.c_str() );
    std::ifstream fin ( fname.c_str() );
    int count;
    fin >> count;
    fin.close();
    if ( count > 0 ){
      output = module[i];
      break;
    }
  }

  return output;
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

  int num_threads(8);

  // int tagging=0;

  double alpha=-1;

  std::size_t sz;

  // while ( (opt = getopt ( argc, argv, "i:o:c:m:RXQ" ) ) != -1 ){
  while ( (opt = getopt ( argc, argv, "i:o:c:m:t:RX" ) ) != -1 ){
    switch ( opt ) {
    case 'i': infile=optarg; break; // infile name
    case 'o': prefix=optarg; break; // output prefix
    case 'c': alpha=std::stod(optarg,&sz); break; // cutoff
    case 'm': path=optarg; break; // path to machine learning model
    case 'R': method+="RandomForest"; break; // Random Forest
    case 'X': method+="XGBoost"; break; // XGBoost
    case 't': num_threads=atoi(optarg); break; // XGBoost
    // case 'Q': tagging=1; break; // Tagging SVs instead of ramoving low quality SVs.
    default: std::cout << "ERROR!!! Check options!!!\n\n" ; usage(); return 0 ; 
    }
  }
  
  if ( infile.size() == 0 ){
    std::cout << "ERROR!!! -i (input file) is required\n"; 
    std::cout << "------------------------------------\n";
    std::cout << "\n";
    usage();
    return 0;
  }

  std::ifstream fin ( infile.c_str() );

  if ( ! fin ){
    std::cout << "ERROR!!! There is no input file: " << infile << "\n";
    std::cout << "-------------------------------------------------\n";
    std::cout << "\n";
    usage();
    return 0;
  }

  fin.close();

  if ( method.size() == 0 ){
    method = "RandomForest";
    //method = "XGBoost";
  }
  else if (method == "RandomForestXGBoost" || method == "XGBoostRandomForest") {
    std::cout << "ERROR!!! -R (Random Forest) and -X (XGBoost) can not be used at the same time\n";
    std::cout << "-----------------------------------------------------------------------------\n";
    std::cout << "\n";
    usage();
    return 0;
  }

  if ( alpha < 0 ){
    alpha = 0.4;
  }

  if ( prefix.size() != 0 ) prefix += ".";

  if ( python_version(prefix) != 3 ) {
    std::cout << "ERROR!!! You need python3.\n";
    return 0;
  }
  python_module_check(prefix,method);



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

  std::string err_fname = prefix + "score.err";

  std::cout << "Calculating features\n";
  calc_feature(infile,feature_file);

  std::cout << "Scoring_command:\t";
  //command = "OMP_NUM_THREADS=4 python3 scorer_" + method + " " + feature_file + " " + score_file + " " + path + " > " + err_fname + " 2>&1" ;
  command = "python scorer_" + method + " " + feature_file + " " + score_file + " " + path + " > " + err_fname + " 2>&1" ;
  echo="echo \"" + command + "\"";

  system ( echo.c_str() );
  system ( command.c_str() );

  std::string check_none;
  fin.open(err_fname.c_str());
  while ( fin >> check_none );
  fin.close();

  if ( check_none != "None" ){
    std::cout << "ERROR!!! Python module of etching_sorter did not run properly\n";
    return 0;
  }

  std::cout << "Removing false positives\n";
  // razor ( infile, score_file, outfile, alpha, method, tagging);
  razor ( infile, score_file, prefix, alpha, method);

  return 0;
}


