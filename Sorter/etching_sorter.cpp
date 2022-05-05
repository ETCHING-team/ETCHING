//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#include "etching_sorter.hpp"
#include <stdlib.h>
#include <random>
#include <sstream>

void usage(){
  std::cout
    << "Usage:\tetching_sorter [options] -i input.vcf -o output_prefix [options]\n"
    << "\n"
    << "Required:\n"
    << "\t" << "-i (string)\t" << "input vcf file\n"
    << "\t" << "-p (string)\t" << "Prefix of machine learning model files\n"
    << "\n"
    << "Options:\n"
    << "\t" << "-o (string)\t" << "Prefix of output vcf file\n"
    << "\t" << "-c (double)\t" << "Cut-off parameter [0.4]\n"
    << "\t" << "-t (int)   \t" << "Number of threads [8]\n"
    << "\t" << "-R         \t" << "Random Forest [default]\n"
    << "\t" << "-X         \t" << "XGBoost\n"
    << "\t" << "NOTE: Machine learning method must be matched with model files.\n"
    << "\n";
}


int make_score_file ( std::string, std::string, std::string, std::string);
  
int main ( int argc , char ** argv ){
  if ( argc == 1 ){
    usage();
    return 1;
  }

  int opt;

  std::string infile;
  std::string prefix;
  std::string method;
  std::string method_R;
  std::string method_X;
  std::string method_L;
  std::string ML_prefix;

  int num_threads(8);

  double cutoff=-1;

  std::size_t sz;

  while ( (opt = getopt ( argc, argv, "i:o:c:t:p:RXLh" ) ) != -1 ){
    switch ( opt ) {
    case 'i': infile=optarg; break; // infile name
    case 'o': prefix=optarg; break; // output prefix
    case 'c': cutoff=std::stod(optarg,&sz); break; // cutoff
    case 't': num_threads=atoi(optarg); break; // number of threads
    case 'p': ML_prefix=optarg; break; // prefix of ML model files
    case 'R': method_R="RandomForest"; break;
    case 'X': method_X="XGBoost"; break;
    case 'L': method_L="LightGBM"; break;
    case 'h': usage() ; return 1;
    default:
      {
	std::cout << "ERROR!!! Check options!!!\n\n";
	std::cout << (char) opt << "\n\n";
	usage();
	return 2 ;
      }
    }
  }
  
  if ( infile.size() == 0 ){
    std::cout << "ERROR!!! -i (input file) is required\n"; 
    std::cout << "------------------------------------\n";
    std::cout << "\n";
    usage();
    return 3;
  }

  std::ifstream fin ( infile.c_str() );

  if ( ! fin.good() ){
    std::cout << "ERROR!!! There is no input file: " << infile << "\n";
    std::cout << "-------------------------------------------------\n";
    std::cout << "\n";
    usage();
    return 4;
  }

  fin.close();

  if ( ML_prefix.size() == 0 ){
    std::cout << "ERROR!!! -m option is required (the prefix of machine learning model)\n";
    std::cout << "---------------------------------------------------------------------\n";
    std::cout << "\n";
    usage();
    return 4;
  }

  if ( method_R.size() == 0 && method_X.size() == 0 && method_L.size() == 0 ){
    method = "RandomForest";
  }
  else if ( method_R.size() != 0 && method_X.size() == 0 && method_L.size() == 0 ){
    method = "RandomForest";
  }
  else if ( method_R.size() == 0 && method_X.size() != 0 && method_L.size() == 0 ){
    method = "XGBoost";
  }
  else if ( method_R.size() == 0 && method_X.size() == 0 && method_L.size() != 0 ){
    method = "LightGBM";
  }
  else {
    std::cout << "ERROR!!! Please use only one of -R (Random Forest), or -X (XGBoost).\n";
    std::cout << "-----------------------------------------------------------------------------\n";
    std::cout << "\n";
    usage();
    return 5;
  }


  if ( cutoff < 0 ){
    cutoff = 0.4;
  }

  if ( prefix.size() != 0 ) prefix += ".";


  std::string outfile = prefix + "etching_sorter.vcf";
  std::string score_file = prefix + "score.txt";

  std::cout << "Method\t" << method << "\n";
  std::cout << "Cut-off\t" << cutoff << "\n";
  std::cout << "Out-file\t" << outfile << "\n";

  std::cout << "Scoring_command:\t";
  if ( method == "RandomForest" ){
    std::cout << "Print feature file\n";
    std::string feature_file = prefix + "feature_table.rf.txt";
    print_rf_feature(infile,feature_file);

    for ( int i = 1 ; i <= 10 ; i ++ ){
      std::cout << "\n";
      std::cout << "---------------------------------\n";
      std::cout << "\n";
      std::string model_file = ML_prefix + "_" + std::to_string(i) + ".forest";
      std::string sub_prefix = prefix + "rf_" + std::to_string(i) ;
      std::string command 
	= "ranger --verbos --treetype 3 --splitrule 4 --nthreads " + std::to_string(num_threads)
	+ " --outprefix " + sub_prefix 
	+ " --file " + feature_file
	+ " --predict " + model_file ;
      std::cout << command << " 2>&1\n";
      int status = system ( command.c_str() );
      if ( status != 0 ){
	std::cout << "\n";
	std::cout << "---------------------------------\n";
	std::cout << "\n";
	std::cout << "ERROR!!! ranger did run properly.\n";
	exit(EXIT_FAILURE);
      }
    }
  }
  else if ( method == "XGBoost" ){
    std::string config_file = prefix + "xgb.conf";
    std::ofstream fout ( config_file.c_str() );
    fout << "task=pred\nbooster=gbtree\nobjective=binary:logistic\neta=0.2\ngamma=3\nmax_depth=30\nnum_round=500\nsave_period=0\n";
    fout.close();
    std::cout << "Print feature file\n";
    std::string feature_file = prefix + "feature_table.xgb.txt";
    print_xgb_feature(infile,feature_file);
    for ( int i = 1 ; i <= 10 ; i ++ ){
      std::cout << "\n";
      std::cout << "---------------------------------\n";
      std::cout << "\n";
      std::string model_file = ML_prefix + "_" + std::to_string(i) + ".model";
      std::string sub_prefix = prefix + "xgb_" + std::to_string(i) + ".prediction";
      std::string command 
	= "xgboost " + config_file + " "
	+ "nthread=" + std::to_string(num_threads) + " "
	+ "test:data=" + feature_file + " "
	+ "model_in=" + model_file + " "
	+ "name_pred=" + sub_prefix;
      std::cout << command << " 2>&1 \n";
      int status = system ( command.c_str() );
      if ( status != 0 ){
	std::cout << "\n";
	std::cout << "---------------------------------\n";
	std::cout << "\n";
	std::cout << "ERROR!!! xgboost did run properly.\n";
	exit(EXIT_FAILURE);
      }
    }
  }


  // TODO :: score.file :: format = ( score \t id ) with a header
  int status;
  status=make_score_file ( infile, score_file, prefix, method);
  if (status != 0 ){
    return status;
  }


  // TODO: error check

  std::cout << "Removing false positives\n";
  razor ( infile, score_file, prefix, cutoff, method);

  return 0;
}


int make_score_file ( std::string infile, std::string score_file, std::string prefix, std::string method){
  std::vector < std::string > id_vec;
  std::vector < double > score_sum_vec;
  std::ifstream fin;
  std::string tmp;
  std::string id;

  std::string prediction_prefix = prefix;
  if ( method == "RandomForest" ){
    prediction_prefix += "rf";
  }
  else if ( method == "XGBoost" ){
    prediction_prefix += "xgb";
  }

  fin.open ( infile.c_str() );
  while ( std::getline ( fin , tmp ) ){
    if (tmp[1]!='#') break;
  }

  while ( std::getline ( fin , tmp ) ){
    std::stringstream ss ( tmp );
    ss >> id >> id >> id;
    id_vec.push_back(id);
  }
  fin.close();

  std::size_t Size=id_vec.size();
  score_sum_vec.resize(Size,0);

  for ( std::size_t i = 0 ; i < 10 ; i ++ ){
    std::vector < double > score_vec;
    std::string prediction = prediction_prefix + "_" + std::to_string(i+1) + ".prediction";
    double score;
    fin.open ( prediction.c_str() );
    if ( method == "RandomForest" ) fin >> tmp ;
    while ( fin >> score ){
      score_vec.push_back(score);
    }
    fin.close();
    std::size_t Size1=score_vec.size();
    if ( Size != Size1 ){
      std::cout << "ERROR!!! Number of SV is not matched with the number of score in etching_sorter." << "\n";
      return 6;
    }
    for ( std::size_t j = 0 ; j < Size ; j ++ ){
      score_sum_vec[j] += score_vec[j]/10;
    }
  }


  std::ofstream fout ( score_file.c_str() );
  fout << "Score\tID\n";
  for ( std::size_t i = 0 ; i < Size ; i ++ ){
    fout << score_sum_vec[i] << "\t" << id_vec[i] << "\n";
  }
  fout.close();

  return 0;
}
