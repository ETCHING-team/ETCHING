//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

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
#include <tuple>

#include <api/BamAlignment.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>

#include "etching_info.hpp"
#include "types.hpp"
#include "caller_functions.hpp"

#ifndef VCF_LINE_CLASS
#define VCF_LINE_CLASS

class VCF_LINE{
private:
  void input_features(std::string feature_str);
  void input_features_2(std::string feature_str);
  
public:
  VCF_LINE(){};
  VCF_LINE(std::string);
  ~VCF_LINE(){};

  void parse(std::string);
  void parse_etching(std::string);

  Feature parse_feature(std::string feature);

  std::string chr1;
  pos_type pos1;
  pos_type dir1;

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
  pos_type pos2;
  pos_type dir2;

  std::string strand;
  

  std::string mate_id; // for BND
  //std::string source_bnd; // Source BND of typed SV
  std::string svtype; // INS, DEL, DUP, INV, BND
  pos_type svlen; // SV length
  std::string scoremethod; // One of Random forest or XGBoost

  // features of first mate
  int cr; // number of Clipped Reads supporing the variation
  int sr; // number of Split Reads supporing the variation
  int pe; // number of PE reads supporing the variation
  int mq; // mapping quailty
  double depdif; // depth difference
  int nxa; // number of XA tags (not using for ML)
  int tcb; // Number of total clipped base
  double entropy; // path entorpy (not using for ML)

  // features of second mate
  int cr2; // number of Clipped Reads supporing the variation
  int sr2; // number of Split Reads supporing the variation
  int pe2; // number of PE reads supporing the variation
  int mq2; // mapping quality
  double depdif2; // depth difference
  int nxa2; // number of XA tags (not using for ML)
  int tcb2; // Number of total clipped base
  double entropy2; // path entropy (not using for ML)

  // global features
  double purity; // tumor purity (not using for ML)
  double seqdep; // sequencing depth (not using for ML)

  pos_type read_number;

  void modify_svtype_info();
  void make_info();

  std::string to_string();
  std::string to_string_short();

  void resize_tool_comp(int Size);

};

VCF_LINE return_mate(VCF_LINE input);

#endif

///////////////////////////////////////////////////////////

#ifndef BP_CLASS
#define BP_CLASS

class BP{
private:
public:
  BP(){
    std::get<0>(this->bp)=0;
    std::get<1>(this->bp)=0;
    std::get<2>(this->bp)=0;
  };
  BP(pos_type a){
    std::get<0>(this->bp)=a;
    std::get<1>(this->bp)=0;
    std::get<2>(this->bp)=0;    
  }
  BP(pos_type a, pos_type b){
    std::get<0>(this->bp)=a;
    std::get<1>(this->bp)=b;
    std::get<2>(this->bp)=0;    
  }
  BP(pos_type a, pos_type b, pos_type c){
    this->bp=std::make_tuple (a,b,c);
  }
  BP(BP_type bp){this->bp=bp;}
  ~BP(){};

  BP_type bp;

  pos_type chr() const {return std::get<0>(this->bp);}
  pos_type pos() const {return std::get<1>(this->bp);}
  pos_type dir() const {return std::get<2>(this->bp);}

  void chr(pos_type a) {std::get<0>(this->bp)=a;}
  void pos(pos_type b) {std::get<1>(this->bp)=b;}
  void dir(pos_type c) {std::get<2>(this->bp)=c;}

  void flip(){std::get<2>(this->bp)*=-1;};

  void insert (BP_type bp){this->bp = bp;}

  BP operator + (pos_type x) const {BP a(*this); std::get<1>(a.bp) += x; return a;}
  BP operator - (pos_type x) const {BP a(*this); std::get<1>(a.bp) -= x; return a;}

  BP & operator += (pos_type x) {std::get<1>(this->bp) += x; return *this;}
  BP & operator -= (pos_type x) {std::get<1>(this->bp) -= x; return *this;}

  BP & operator =  (BP_type bp) {this->bp = bp;return *this;}
  
  bool operator == (BP input) const {return this->bp==input.bp;}
  bool operator != (BP input) const {return this->bp!=input.bp;}
  bool operator <  (BP input) const {return this->bp< input.bp;}
  bool operator >  (BP input) const {return this->bp> input.bp;}
  bool operator <= (BP input) const {return this->bp<=input.bp;}
  bool operator >= (BP input) const {return this->bp>=input.bp;}
};


#endif

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef BND_CLASS
#define BND_CLASS

class BND{
private:
public:
  BND(){};
  BND(BP first, BP second){
    this->first = first;
    this->second = second;
  };
  BND(BND_type i) {this->first=i.first;this->second=i.second;};
  ~BND(){};

  BP first;
  BP second;

  const bool operator == (BND input) const {return (first.bp==input.first.bp&&second.bp==input.second.bp);}
  const bool operator != (BND input) const {return ! (*this == input) ;}
  const bool operator <  (BND input) const {
    return this->first<input.first?1:this->first>input.first?0:this->second<input.second;
  }
  const bool operator >  (BND input) const {
    return this->first>input.first?1:this->first<input.first?0:this->second>input.second;
  }
  const bool operator <= (BND input) const {return ! (*this < input);}
  const bool operator >= (BND input) const {return ! (*this > input);}

  Position Pos1() const {
    return std::make_pair(this->first.chr(),this->first.pos());
  }

  Position Pos2() const {
    return std::make_pair(this->second.chr(),this->second.pos());
  }

  pos_type chr1() const {return this->first.chr();}
  pos_type pos1() const {return this->first.pos();}
  pos_type dir1() const {return this->first.dir();}
  pos_type chr2() const {return this->second.chr();}
  pos_type pos2() const {return this->second.pos();}
  pos_type dir2() const {return this->second.dir();}

};

#endif

////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef VCF_CLASS
#define VCF_CLASS

// typedef std::map < BP , std::map < BP , VCF_LINE > > BP_MAP;
typedef std::map < BND , VCF_LINE > VCF_MAP;

class VCF{
private:
  std::vector < std::pair < std::string , std::size_t > > genome_info;
  void build_id_ref_map(const std::string infile);
  char GetNucl(const std::string & , const std::size_t &);
  std::string get_genome_sequence_id ( std::string input );
  
public:
  VCF();
  VCF(const std::string input_file);
  ~VCF();

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
  pos_type insert_size;
  pos_type read_length;

  BamTools::RefVector refvector;
  std::map < std::string , int > id_ref_map;
  std::map < int , std::string > ref_id_map;
  std::map < std::string , std::string > genome;

  VCF_MAP::iterator begin();
  VCF_MAP::iterator end();
  VCF_MAP::iterator find (BND bnd) ;
  VCF_MAP::iterator find (BP bp1, BP bp2) ;

  std::size_t size();
  std::map<Position,double> calc_link_entropy (const std::string pair_file);

  void read_vcf_file(const std::string infile);
  void get_genome();
  void get_genome(std::string);
  void make_header();
  void make_header_short();
  void clear();
  void write();
  void write_no_header();
  void fwrite(std::string);
  void write_short();
  void fwrite_short(std::string);

  void insert ( std::string chr1, pos_type pos1, std::string chr2, pos_type pos2, std::string sv_id, std::string mate_id, std::string strand, pos_type sr_val, std::string svtype);
  void insert ( VCF_LINE vcf_line );

  void calc_features(const std::string input_bam, const pos_type read_length, const pos_type insert_size, const pos_type confi_window, const pos_type pe_window_size);

  void add_features_in_id();

  bool check_vcf( BND bnd ) ;
  bool check_vcf( BP bp1, BP bp2 ) ;
  bool check_vcf( std::string chr1, pos_type pos1, std::string chr2, pos_type pos2, std::string strand ) ;
  bool check_vcf( pos_type chr1, pos_type pos1, pos_type chr2, pos_type pos2, std::string strand ) ;
  // bool check_vcf_bidir( BND bnd ) ;
  bool check_vcf_bidir( BP bp1, BP bp2 ) ;
  bool check_vcf_bidir( std::string chr1, pos_type pos1, std::string chr2, pos_type pos2 ) ;
  bool check_vcf_bidir( pos_type chr1, pos_type pos1, pos_type chr2, pos_type pos2 ) ;

  void make_info();
  void get_reference_information ( std::vector < std::string > file_vec );

  VCF_LINE & operator [](BND bnd);
};

double return_depdif(Position Pos, std::vector < double > & dep_vec, pos_type Size, double read_length);
void copy_info ( VCF & source , VCF & target );
VCF typing_SV(VCF & source);
void calc_features(const std::string input_bam, std::vector < VCF > & container_vec, const pos_type read_length, const pos_type insert_size, const pos_type confi_window) ;

#endif

