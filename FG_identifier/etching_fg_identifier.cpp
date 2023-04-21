//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#include <string>
#include <iostream>
#include <unistd.h>
#include <getopt.h>

#include "etching_info.hpp"
#include "find_fg.hpp"
#include "make_bp_pair.hpp"
#include "make_gene_list.hpp"

void fg_identifier_usage(int argc, char ** argv){
  std::cout
    << "Program: " << get_program_name(argc,argv) << "\n"
    << "Version: " << ETCHING_VERSION << " (" << RELEASE_DATE << ")\n"
    << "Usage: etching_fg_identifier [options] -i input.vcf -a annotation.gtf [options]\n"
    << "\n"
    << "[Options]\n"
    << "-i (string) \tinput.vcf (required)\n"
    << "-o (string) \toutput.txt (required)\n"
    << "-a (string) \tannotation.gtf (required)\n"
    << "-w (int)    \tWindow size for fusion-gene prediction [10000]\n"
    << "-f (string) \tField name [gene_name]\n"
    << "-s          \tPredict FGs considering strand\n"
    << "\n"
    << "-h          \tPrint this message\n";
}

int main ( int argc , char ** argv ){

  if ( argc == 1 ){
    fg_identifier_usage(argc, argv);
    return 0;
  }

  std::string input_vcf_file;
  std::string outfile;
  std::string input_annotation_file;
  int window=10000;
  std::string fieldName = "gene_name";
  int strand_aware=0;
  int prom_win=0;         // for debugging
  int opt = 0;

  while (1){
    static struct option long_options[]=
      {
        {"strand-aware",     no_argument,       &strand_aware,     's'},
        {"input-vcf",        required_argument, 0,                 'i'},
        {"output-file",      required_argument, 0,                 'o'},
        {"annotation-gtf",   required_argument, 0,                 'a'},
        {"field-name",       required_argument, 0,                 'f'},
        {"fusion-window",    required_argument, 0,                 'w'},
        {"help",             no_argument,       0,                 'h'},
        {0,                  0,                 0,                 0}
      };

    int opt_ind=0;
    opt=getopt_long(argc,argv,"sri:o:a:f:w:h",long_options,&opt_ind);
    if (opt==-1) break;

    switch(opt){
    case 's': strand_aware = 1; break;
    case 'i': input_vcf_file = optarg; break;
    case 'o': outfile = optarg; break;
    case 'a': input_annotation_file = optarg; break;
    case 'w': window = atoi(optarg); break;
    case 'f': fieldName = optarg; break;
    case 'h': fg_identifier_usage(argc,argv) ; return 0;
    default: std::cout << "ERROR!!! Invalid option: " << optarg << "\n" ; return -1;
    }
  }

  if ( input_vcf_file.size() == 0 ){
    std::cout << "ERRORR!!! -i option is required.\n"
	      << "-----------------------------------\n";
    fg_identifier_usage(argc,argv);
    return -1;
  }

  if ( input_annotation_file.size() == 0 ){
    std::cout << "ERRORR!!! -a option is required.\n"
	      << "-----------------------------------\n";
    fg_identifier_usage(argc,argv);
    return -1;
  }

  std::string bp_pair = input_vcf_file + ".bp_pair" ;
  std::string gene_list = input_annotation_file + ".gene_list" ;

  std::string SW[2];
  SW[0]="NO";
  SW[1]="YES";

  std::cerr << "[Input file: " << input_vcf_file << "]\n";
  std::cerr << "[Annotation: " << input_annotation_file << "]\n";
  std::cerr << "[Field name: " << fieldName<< "]\n";
  std::cerr << "[Window size: " << window << "]\n";
  std::cerr << "[Strand-ware: " << SW[strand_aware] << "]\n";

  std::ifstream fin( input_vcf_file.c_str() );
  if ( fin.good() ){
    make_bp_pair(input_vcf_file,bp_pair);
  }
  else {
    std::cout << "ERROR!!! Please check input file: -i " << input_vcf_file << "\n";
    return -1;
  }
  fin.close();

  fin.open( input_annotation_file.c_str() );
  if ( ! fin.good() ){
    std::cout << "ERROR!!! Please check input file: -a " << input_annotation_file << "\n";
    return 1;
  }
  fin.close();


  fin.open( gene_list.c_str() );
  if ( ! fin.good() ){
    make_gene_list(input_annotation_file, "gene", fieldName, gene_list, prom_win);
  }
  fin.close();

  std::stringstream fg_output = find_fg ( bp_pair, gene_list, window, strand_aware );
  std::ofstream fout ( outfile.c_str() );
  if ( ! fout.good() ){
    std::cerr<<"ERROR!!!"<<__FILE__<<"::"<<__FUNCTION__<<"::"<<__LINE__<<"::Out file cannot be open::"<<outfile<<"\n";
    exit (EXIT_FAILURE);
  }
  fout << fg_output.str();
  fout.close();

  return 0;
}




