//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#ifndef ETCHING_SORTER
#define ETCHING_SORTER

#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <map>

#define NO_OF_FEATURE 6


using Feature = std::map < std::string , double > ;
using FEATURE = std::map < std::string , Feature > ;


Feature harmonic_average (Feature f1, Feature f2){
  for ( auto & i : f1 ){
    if ( i.second * f2[i.first] != 0 ){
      i.second = ( 2 * i.second * f2[i.first] ) / ( i.second + f2[i.first]) ;
    }
    else i.second = 0;
  }
  return f1;
}


std::string cut_str(std::string & input, std::string key){
  std::size_t found = input.find(key);
  std::string output;
  if ( found != std::string::npos ){
    output = input.substr(0,found);
    input = input.substr(found+key.size());
  }
  else{
    output = input;
    input.clear();
  }
  return output;
}                                    

std::string return_svtype(std::string info){
  std::string key = "SVTYPE=";
  std::size_t found = info.find(key);
  return info.substr(found+key.size(),3);
} 

Feature parse_feature ( std::string format , std::string format1){
  Feature output;

  std::string name;
  std::string tmp;
  double value;
  std::string key=":";
  std::size_t sz;

  while ( format.size() != 0 ){
    name = cut_str(format,key);
    tmp = cut_str(format1,key);
    value = std::stod(tmp,&sz);
    output[name]=value;
  }

  return output;
}



void print_header(std::ofstream & fout, std::vector < std::string > feat_vec){
  for ( int i = 0 ; i < NO_OF_FEATURE ; i ++ ){
    fout << feat_vec[i] << "\t";
  }
  fout << "TorF\n";
}



void print_rf_feature ( std::string infile, std::string outfile ){

  std::string tmp;
  std::string id;
  std::string info;
  std::string format;
  std::string format1;
  std::string format2;

  Feature feat;
  Feature feat1;
  Feature feat2;

  std::vector < std::string > feat_vec(NO_OF_FEATURE);
  feat_vec[0]="CR";
  feat_vec[1]="SR";
  feat_vec[2]="PE";
  feat_vec[3]="MQ";
  feat_vec[4]="DEPDIF";
  feat_vec[5]="TCB";

  std::ifstream fin ( infile );
  std::ofstream fout ( outfile );
  while ( std::getline ( fin , tmp ) ){
    if ( tmp[1] != '#' ){
      break;
    }
  }

  print_header(fout, feat_vec);

  while ( fin >> tmp >> tmp >> id >> tmp >> tmp >> tmp >> tmp >> info >> format >> format1 >> format2 ){

    feat1 = parse_feature(format,format1);
    feat2 = parse_feature(format,format2);

    feat = harmonic_average ( feat1 , feat2 );

    for ( int i = 0 ; i < NO_OF_FEATURE ; i ++ ){
      fout << feat[feat_vec[i]] << "\t" ;
    }
    
    fout << "1\n";
  }
  fin.close();
  fout.close();
}


#include <iomanip>

void print_xgb_feature ( std::string infile, std::string outfile ){

  std::string tmp;
  std::string id;
  std::string info;
  std::string format;
  std::string format1;
  std::string format2;

  Feature feat;
  Feature feat1;
  Feature feat2;

  std::vector < std::string > feat_vec(NO_OF_FEATURE);
  feat_vec[0]="CR";
  feat_vec[1]="SR";
  feat_vec[2]="PE";
  feat_vec[3]="MQ";
  feat_vec[4]="DEPDIF";
  feat_vec[5]="TCB";

  std::ifstream fin ( infile );
  std::ofstream fout ( outfile );
  while ( std::getline ( fin , tmp ) ){
    if ( tmp[1] != '#' ){
      break;
    }
  }

  while ( fin >> tmp >> tmp >> id >> tmp >> tmp >> tmp >> tmp >> info >> format >> format1 >> format2 ){


    feat1 = parse_feature(format,format1);
    feat2 = parse_feature(format,format2);

    feat = harmonic_average ( feat1 , feat2 );
    fout.precision(6);
    fout << "1 ";
    for ( int i = 0 ; i < NO_OF_FEATURE - 1 ; i ++ ){
      // fout << std::fixed << std::setprecision(8) << i << ":" << (double) feat[feat_vec[i]] << "\t" ;
      fout << i << ":" << (double) feat[feat_vec[i]] << "\t" ;
    }
    fout << NO_OF_FEATURE - 1 << ":" << (double) feat[feat_vec[NO_OF_FEATURE-1]] << "\n" ;

  }
  fin.close();
  fout.close();
}



void razor ( std::string infile, std::string score_file, std::string prefix, double cutoff, std::string method){

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

  std::string outfile = prefix + "etching_sorter.vcf";
  std::ofstream fout ( outfile.c_str() );
  
  std::string scoredfile = prefix + "unfiltered.vcf";
  std::ofstream fout_score ( scoredfile.c_str() );

  fin.open(infile);
  while ( std::getline(fin,tmp) ){
    fout << tmp <<"\n";
    fout_score << tmp <<"\n";
    if ( tmp[1] != '#' ) break;
  }

  while ( fin >> vcf_content[0] ){
    for ( int i = 1 ; i < 11 ; i ++ ){
      fin >> vcf_content[i] ;
    }
    if ( score_map[vcf_content[2]] >= cutoff ){
      vcf_content[6] = "PASS";
    }
    else {
      vcf_content[6] = "LOWQUAL";
    }
    vcf_content[5] = std::to_string(score_map[vcf_content[2]]);
    vcf_content[7] += ";SVSCORE=" + std::to_string(score_map[vcf_content[2]]) + ";SCOREMETHOD=" + method;
    // vcf_content[7] += ";SCOREMETHOD=" + method;
    for ( int i = 0 ; i < 10 ; i ++ ){
      fout_score << vcf_content[i] << "\t";
    }
    fout_score << vcf_content[10] << "\n";
    if ( vcf_content[6] == "PASS" ){
      for ( int i = 0 ; i < 10 ; i ++ ){
	fout << vcf_content[i] << "\t";
      }
      fout << vcf_content[10] << "\n";
    }
  }

  fin.close();
  fout.close();
  fout_score.close();
}

#endif
