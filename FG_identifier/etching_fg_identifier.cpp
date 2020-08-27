//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#include <string>
#include <iostream>
#include "find_fg.hpp"
#include "make_bp_pair.hpp"
#include "make_gene_list.hpp"

int main ( int argc , char ** argv ){
  if ( argc == 1 ){
    std::cout << "Usage:  etching_fg_identifier  input.vcf  annotation.gtf  [field_name (default: gene_name)]\n";
    return 0;
  }
  
  std::string bp_pair = argv[1] ;
  bp_pair += ".bp_pair";
  std::string gene_list = argv[2] ;
  gene_list += ".gene_list";
  std::string fieldName = "gene_name";
  if ( argc == 4 ) fieldName = argv[3];

  std::ifstream fin( argv[1] );
  if ( fin.good() ){
    make_bp_pair(argv[1],bp_pair);
  }
  else {
    std::cout << "ERROR!!! Please check input file: " << argv[1] << "\n";
    return 0;
  }
  fin.close();

  fin.open( argv[2] );
  if ( ! fin.good() ){
    std::cout << "ERROR!!! Please check input file: " << argv[2] << "\n";
    return 0;
  }
  fin.close();


  fin.open( gene_list.c_str() );
  if ( ! fin )
    make_gene_list(argv[2], "gene", fieldName ,gene_list);

  fin.close();

  std::cout << "#gene_1\tgene_2\tchr_1\tpos_1\tstrand_1\tchr_2\tpos_2\tstrand_2\n";

  find_fg(bp_pair,gene_list);

  return 0;
}
