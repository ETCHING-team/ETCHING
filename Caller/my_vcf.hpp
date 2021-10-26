//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------


#ifndef MY_VCF
#define MY_VCF

#include <string>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <assert.h>
#include <iterator>
#include <time.h>
#include <cmath>

#include <api/BamAlignment.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>

#include "types.hpp"
#include "caller_functions.hpp"



class VCF{
private:
  void input_features(std::string feature_str);
  void input_features_2(std::string feature_str);
  std::string cut_str(std::string & input, std::string key);
  
public:
  VCF(){};
  VCF(std::string);
  ~VCF(){};

  std::string parse(std::string);
  void parse_etching(std::string);

  Feature parse_feature(std::string feature);

  std::string chr1;
  int pos1;
  std::string sv_id;
  std::string sv_id_add;
  std::string ref;
  std::string alt;
  std::string qual;
  std::string filter;
  std::string info;
  
  std::string feature_str;
  std::string feature_str1;
  std::string feature_str2;

  Feature feature;
  Feature feature1;
  Feature feature2;

   // for comparing with other tools
  std::vector < int > tool_comp;
  int tool_count;
  
  std::string chr2;
  int pos2;

  std::string strand;


  std::string mate_id; // for BND
  std::string svtype; // INS, DEL, DUP, INV, BND
  int svlen; // SV length
  std::string scoremethod;

  // features of first mate
  int cr; // number of Clipped Reads supporing the variation
  int sr; // number of Split Reads supporing the variation
  int pe; // number of PE reads supporing the variation
  int mq;
  double depdif;
  int nxa;
  int tcb; // Number of total clipped base
  double entropy;

  // features of second mate
  int cr2; // number of Clipped Reads supporing the variation
  int sr2; // number of Split Reads supporing the variation
  int pe2; // number of PE reads supporing the variation
  int mq2;
  double depdif2;
  int nxa2;
  int tcb2; // Number of total clipped base
  double entropy2;

  // global features
  double purity;
  double seqdep;

  int read_number;


  void modify_svtype_info();
  void make_info();

  std::string to_string();
  std::string to_string_short();

  void resize_tool_comp(int Size);

};

VCF return_mate(VCF input);


///////////////////////////////////////////////////////////
using VCF_MAP=std::map < Position , std::vector < VCF > >;


class VCF_CLASS{
private:
  std::vector < std::pair < std::string , std::size_t > > genome_info;
  void build_id_ref_map(const std::string infile);
  char GetNucl(const std::string & , const std::size_t &);
  std::string get_genome_sequence_id ( std::string input );

public:
  VCF_CLASS();
  VCF_CLASS(const std::string input_file);
  ~VCF_CLASS();

  std::string etching_version="ETCHING_v1.3.7 (2021.10.26.)";
  
  // Main container
  VCF_MAP vcf_map;
  
  std::string vcf_file;
  std::string reference;

  std::string metainfo;
  std::string header;

  double seqdep;
  double purity;
  std::string seqtype;

  std::string bam_file;
  std::string single_file;
  std::string pair_file;
  int insert_size;
  int read_length;

  
  BamTools::RefVector refvector;
  std::map < std::string , int > id_ref_map;
  std::map < int , std::string > ref_id_map;

  std::map < std::string , std::string > genome;

  void read_vcf_file(const std::string infile);

  void get_genome();
  void get_genome(std::string);
  void make_header();
  void make_header_short();
  void clear();

  VCF_MAP::iterator begin();
  VCF_MAP::iterator end();

  void write();
  void write_no_header();
  void fwrite(std::string);
  void write_short();
  void fwrite_short(std::string);
  std::size_t size();

  void insert ( std::string chr1, int pos1, std::string chr2, int pos2, std::string sv_id, std::string mate_id, std::string strand, int sr_val, std::string svtype);


  void insert ( VCF vcf );

  void calc_features(const std::string input_bam, const int read_length, const int insert_size, const int confi_window);
  
  void add_features_in_id();

  VCF_MAP::iterator find(Position Pos);
  bool check_vcf( std::string chr1, int pos1, std::string chr2, int pos2, std::string BND_type );
  std::map<Position,double> calc_link_entropy (const std::string pair_file);

  std::vector < VCF > & operator [](Position Pos);

  void make_info();
};
  
double return_depdif(Position Pos, std::vector < double > & dep_vec, int Size, double read_length);

void copy_info ( VCF_CLASS & source , VCF_CLASS & target );

VCF_CLASS typing_SV(VCF_CLASS & source);

void calc_features(const std::string input_bam, std::vector < VCF_CLASS > & container_vec,
		   const int read_length, const int insert_size, const int confi_window) ;

#endif

