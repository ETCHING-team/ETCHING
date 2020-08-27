//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <limits>
#include <cmath>

int main ( int argc , char ** argv ){
  if ( argc != 4 ){
    std::cout << "estimate_coverage  TUMOR_DB_name  read_length  k-mer_length\n";
    return 0;
  }

  std::string tumor = argv[1];
  double read_length = (double) atoi ( argv[2] );
  double kl = (double) atoi ( argv[3] );


  std::vector <std::size_t> depth;
  std::vector <std::size_t> freq;
  std::size_t d,f;
  // std::size_t local_min;
  std::size_t local_max;

  double coverage;

  std::string command = "kmc_tools transform " + tumor + " histogram " + tumor + ".hist";
  system(command.c_str());

  std::string tumor_hist = tumor + ".hist";

  std::ifstream fin(tumor_hist.c_str());
  while ( fin >> d >> f ){
    depth.push_back(d);
    freq.push_back(f);
  }
  fin.close();
  

  {
    size_t i = 1 ; 
    for ( ; i < depth.size() ; i ++ ){
      if ( freq[i-1] < freq[i] ){
	// local_min = depth[i-1];
	break;
      }
    }
    for ( ; i < depth.size() ; i ++ ){
      if ( freq[i-1] > freq[i] ){
	local_max = depth[i-1];
	break;
      }
    }
  }

  coverage = (double) local_max * read_length;
  coverage /= ( (double) read_length - kl + 1 );
  
  std::cout << std::to_string((int)round(coverage)) << "\n";

  return 0;
}
