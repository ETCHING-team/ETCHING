//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#include <string>
#include <iostream>
#include <unistd.h>

#include "find_fg.hpp"
#include "make_bp_pair.hpp"
#include "make_gene_list.hpp"

void fg_identifier_usage(){
  std::cout << "Usage: etching_fg_identifier [options] -i input.vcf -a annotation.gtf [options]\n"
	    << "\n"
	    << "[Options]\n"
	    << "-f (string)\tField name [gene_name]\n"
	    << "-w (int)   \tWindow [10000]\n"
	    << "-h         \tPrint this message\n";
}

int main ( int argc , char ** argv ){

  if ( argc == 1 ){
    fg_identifier_usage();
    return 0;
  }

  std::string infile;
  std::string annotation;
  int window=10000;
  std::string fieldName = "gene_name";

  int opt = 0;
  while ((opt = getopt ( argc , argv , "i:a:f:w:h") ) != -1 ){
    switch(opt){
    case 'i': infile = optarg; break;
    case 'a': annotation = optarg; break;
    case 'w': window = atoi(optarg); break;
    case 'f': fieldName = optarg; break;
    case 'h': fg_identifier_usage() ; return 0;
    default: std::cout << "ERROR!!! Invalid option: " << optarg << "\n" ; return -1;
    }
  }

  if ( infile.size() == 0 ){
    std::cout << "ERRORR!!! -i option is required.\n"
	      << "-----------------------------------\n";
    fg_identifier_usage();
    return -1;
  }

  if ( annotation.size() == 0 ){
    std::cout << "ERRORR!!! -a option is required.\n"
	      << "-----------------------------------\n";
    fg_identifier_usage();
    return -1;
  }

  std::string bp_pair = infile + ".bp_pair" ;
  // bp_pair += ".bp_pair";
  std::string gene_list = annotation + ".gene_list" ;
  // gene_list += ".gene_list";

  

  std::ifstream fin( infile.c_str() );
  if ( fin.good() ){
    make_bp_pair(infile,bp_pair);
  }
  else {
    std::cout << "ERROR!!! Please check input file: " << infile << "\n";
    return 0;
  }
  fin.close();

  fin.open( annotation.c_str() );
  if ( ! fin.good() ){
    std::cout << "ERROR!!! Please check input file: " << annotation << "\n";
    return 0;
  }
  fin.close();


  fin.open( gene_list.c_str() );
  if ( ! fin )
    make_gene_list(annotation, "gene", fieldName ,gene_list);

  fin.close();

  std::cout << "#gene_1\tgene_2\tchr_1\tpos_1\tstrand_1\tchr_2\tpos_2\tstrand_2\n";

  find_fg ( bp_pair, gene_list, window);

  return 0;
}
