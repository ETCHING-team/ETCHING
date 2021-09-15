//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#include <string>
#include <iostream>
#include <unistd.h>
#include <getopt.h>

#include "find_fg.hpp"
#include "make_bp_pair.hpp"
#include "make_gene_list.hpp"

void fg_identifier_usage(){
  std::cout << "Usage: etching_fg_identifier [options] -i input.vcf -a annotation.gtf [options]\n"
	    << "\n"
	    << "[Options]\n"
	    << "-f (string)\tField name [gene_name]\n"
	    << "-w (int)   \tWindow size for fusion-gene prediction [10000]\n"
	    << "--fusion-window (int)\n"
	    << "           \tWindow size to detect fusion-gene [10000]\n"
	    << "--strand-aware\n"
	    << "           \tPredict FGs connecting nearest genes to a BP-pair awaring strand\n"
            << "           \t[Default: NONE (all possible gene-pairs both sides of each BP-pair)]\n\n"
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
  int strand_aware=0;

  int opt = 0;

  while (1){
    static struct option long_options[]=
      {
        {"strand-aware",   no_argument,       &strand_aware, 0},
        {"input-vcf",      required_argument, 0, 'i'},
        {"annotation-gtf", required_argument, 0, 'a'},
        {"field-name",     required_argument, 0, 'f'},
        {"fusion-window",  required_argument, 0, 'w'},
        {"help",           no_argument,       0, 'h'},
        {0, 0, 0, 0}
      };

    int opt_ind=0;
    opt=getopt_long(argc,argv,"i:a:f:w:h",long_options,&opt_ind);
    if (opt==-1) break;

    switch(opt){
    case   0: strand_aware = 1; break;
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
  std::string gene_list = annotation + ".gene_list" ;

  std::string SW[2];
  SW[0]="NO";
  SW[1]="YES";

  std::cerr << "[Input file: " << infile << "]\n";
  std::cerr << "[Annotation: " << annotation << "]\n";
  std::cerr << "[Field name: " << fieldName<< "]\n";
  std::cerr << "[Window size: " << window << "]\n";
  std::cerr << "[Strand-ware: " << SW[strand_aware] << "]\n";

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

  std::cout << "#gene_1\tgene_2\tchr_1\tpos_1\tstrand_1\tchr_2\tpos_2\tstrand_2\tclass\n";

  find_fg ( bp_pair, gene_list, window, strand_aware );

  return 0;
}
