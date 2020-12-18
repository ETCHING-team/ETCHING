//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#ifndef CALC_FEATURE
#define CALC_FEATURE


#include "etching_sorter.hpp"



void print_header(std::ofstream & fout, std::vector < std::string > feat_vec){
  fout << "ID" << "\t" << "SVTYPE" ;
  for ( int i = 0 ; i < 6 ; i ++ ){
    fout << "\t" << feat_vec[i];
  }
  fout << "\n";
}



void print_feat(std::string id, std::string type, Feature f, std::vector < std::string > feat_vec, std::ofstream & fout){
  fout << id << "\t" << type ;

  for ( int i = 0 ; i < 6 ; i ++ ){
    fout << "\t" << f[feat_vec[i]] ;
  }

  fout << "\n";
}




void calc_feature ( std::string infile, std::string outfile ){

  std::string tmp;
  std::string id;
  std::string type;
  std::string info;
  std::string format;
  std::string format1;
  std::string format2;

  Feature feat;
  Feature feat1;
  Feature feat2;

  std::vector < std::string > feat_vec(6);
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

    type=return_svtype(info);

    feat1 = parse_feature(format,format1);
    feat2 = parse_feature(format,format2);

    feat = harmonic_average ( feat1 , feat2 );

    feat["CR"] /= feat["PURITY"] * feat["SEQDEP"];
    feat["SR"] /= feat["PURITY"] * feat["SEQDEP"];
    feat["PE"] /= feat["PURITY"] * feat["SEQDEP"];
    feat["DEPDIF"] /= feat["PURITY"] * feat["SEQDEP"];
    feat["TCB"] /= feat["PURITY"] * feat["SEQDEP"];

    print_feat(id, type, feat, feat_vec, fout);

  }
  fin.close();
  fout.close();
}



#endif
