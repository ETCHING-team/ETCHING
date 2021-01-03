//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#include "my_vcf.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <vector>

void typer_usage(){
  std::cout << "Usage:  etching_typer  input.vcf\n" ;
}


int main ( int argc , char ** argv ){
  if ( argc == 1 ){
  typer_usage();
  return 0;
}

  std::string infile = argv[1];
  // std::string index = argv[2];
  std::string outfile = infile.substr(0,infile.size()-3) + "etching_typer.vcf";

  std::string id;
  // int64_t length;

  std::string tmp;
  std::vector < std::string > metainfo;
  VCF vcf;
  VCF_CLASS container;
  VCF_CLASS container_SV;

  // int64_t count = 0 ;
  
  std::ifstream fin;
  // fin.open( index.c_str() );
  // while ( fin >> id ){
  //   fin >> length;
  //   container.id_ref_map[id] = count;
  //   container.ref_id_map[count] = id;
  //   count ++;
  // }

  // fin.close();
  fin.open( infile.c_str() );

  /////////////////////////////////////////////////////////////////
  //
  // Reading metainfo
  //

  while ( std::getline ( fin , tmp ) ){
    if ( tmp.substr(0,2) == "##" ){
      metainfo.push_back(tmp);
    }
    else break;
  }
  metainfo.push_back(tmp);

  /////////////////////////////////////////////////////////////////
  //
  // Make id number maps
  //

  const std::string key = "##contig=<ID=";
  const std::string len = ",length=";
  std::string Chr;

  
  /////////////////////////////////////////////////////////////////
  //
  // fill VCF container
  //


  while ( std::getline ( fin , tmp ) ){
    vcf.parse_etching(tmp);
    container.insert(vcf);
  }

  fin.close();


  /////////////////////////////////////////////////////////////////
  //
  // Add metainfo
  //

  // found = 0 ;
  for ( std::size_t i = 0 ; i < metainfo.size() - 2 ; i ++ ){
    container.metainfo += metainfo[i] + "\n";
  }

  container.metainfo += metainfo[metainfo.size()-2];
  container.header = metainfo[metainfo.size()-1];

  
  /////////////////////////////////////////////////////////////////
  //
  // Typing SVs
  //
  
  container_SV = typing_SV_general ( container );
  
  /////////////////////////////////////////////////////////////////
  //
  // print result
  //
  
  container_SV.fwrite(outfile);
  return 0;
}
