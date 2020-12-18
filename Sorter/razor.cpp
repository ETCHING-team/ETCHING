//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#ifndef ETCHING_RAZOR
#define ETCHING_RAZOR

#include "etching_sorter.hpp"
#include <map>
#include <vector>

void razor ( std::string infile, std::string score_file, std::string outfile, double cutoff, std::string method, int tagging){

  std::string id;
  double value;
  std::string tmp;

  std::map < std::string , double > score_map;
  std::vector < std::string > vcf_content(11);

  std::ifstream fin(score_file.c_str());
  std::getline(fin,tmp);
  while ( fin >> value >> id ){
    score_map[id]=value;
  }
  fin.close();

  std::ofstream fout ( outfile.c_str() );
  fin.open(infile);
  while ( std::getline(fin,tmp) ){
    fout << tmp <<"\n";
    if ( tmp[1] != '#' ) break;
  }
  if ( tagging == 0 ){
    while ( fin >> vcf_content[0] ){
      for ( int i = 1 ; i < 11 ; i ++ ){
	fin >> vcf_content[i] ;
      }
      if ( score_map[vcf_content[2]] >= cutoff ){
	vcf_content[6] = "PASS";
	vcf_content[7] += ";SVSCORE=" + std::to_string(score_map[vcf_content[2]]) + ";SCOREMETHOD=" + method;
	for ( int i = 0 ; i < 10 ; i ++ ){
	  fout << vcf_content[i] << "\t";
	}
	fout << vcf_content[10] << "\n";
      }
    }
  }
  else{
    while ( fin >> vcf_content[0] ){
      for ( int i = 1 ; i < 11 ; i ++ ){
	fin >> vcf_content[i] ;
      }
      if ( score_map[vcf_content[2]] >= cutoff ){
	vcf_content[6] = "PASS";
      }
      else {
	vcf_content[6] = "LowQual";
      }
      vcf_content[7] += ";SVSCORE=" + std::to_string(score_map[vcf_content[2]]) + ";SCOREMETHOD=" + method;
      for ( int i = 0 ; i < 10 ; i ++ ){
	fout << vcf_content[i] << "\t";
      }
      fout << vcf_content[10] << "\n";
    }
  }
  fin.close();
  fout.close();
}

#endif
