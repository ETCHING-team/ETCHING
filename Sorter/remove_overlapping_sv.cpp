//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------


#ifndef ETCHING_SOMATIC_FILTER
#define ETCHING_SOMATIC_FILTER

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
#include <sys/time.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "CPU_TIME_measure.hpp"
#include "Peak_memory_measure.hpp"
#include "etching_info.hpp"

using Position = std::pair < int , int > ;

class vcf_record{
private:
  //void input_features(std::string feature_str);
  //void input_features_2(std::string feature_str);
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

  
public:
  vcf_record(){};
  vcf_record(std::string line) { parse(line); }
  ~vcf_record(){};

  std::string chr1;
  int pos1;
  std::string sv_id;
  std::string sv_id_add;
  std::string ref;
  std::string alt;
  std::string qual;
  std::string filter;
  std::string info;
  
  std::string chr2;
  int pos2;

  std::string strand;


  std::string parse(std::string line){

    std::string tmp;
    std::string key;
    std::size_t sz;
    std::size_t found;
  
    chr2 = "";
    pos2 = -1;


    // chr
    found = line.find('\t');
    chr1 = line.substr(0,found);
    line = line.substr(found+1);
  
    // pos
    found = line.find('\t');
    pos1 =  atoi(line.substr(0,found).c_str());
    line = line.substr(found+1);

    // id
    found = line.find('\t');
    sv_id =  line.substr(0,found);
    line = line.substr(found+1);

    // reference sequence
    found = line.find('\t');
    ref =  line.substr(0,found);
    line = line.substr(found+1);

    // alt sequence
    found = line.find('\t');
    alt =  line.substr(0,found);
    line = line.substr(found+1);

    // quality
    found = line.find('\t');
    if ( isdigit(line.substr(0,found)[0] ) ){
      // qual =  atoi(line.substr(0,found).c_str());
      qual =  line.substr(0,found).c_str();
    }
    else{
      qual = ".";
    }
    line = line.substr(found+1);

    // filter
    found = line.find('\t');
    filter =  line.substr(0,found);
    line = line.substr(found+1);

    // information
    info =  line;

    /////////////////////////////////////////////////////////////////
    //
    // read chr2 and pos2
    //

    std::string chr2_key="CHR2=";
    std::string pos2_key=";END=";

    tmp = alt;
    std::string tmp1;

    key = "[";
    tmp1 = cut_str(tmp, key);
    if ( tmp.size() > 0 ){
      chr2 = cut_str(tmp,":");
      tmp = cut_str(tmp,"[");
      pos2 = std::stol(tmp);
    }

    tmp = alt;
    key = "]";
    tmp1 = cut_str(tmp, key);
    if ( tmp.size() > 0 ){ 
      chr2 = cut_str(tmp, ":");
      tmp = cut_str(tmp,key);
      pos2 = std::stol(tmp);
    }


    if ( chr2 == "" ){
      tmp = info;
      tmp1 = cut_str(tmp,chr2_key);
      if ( tmp.size() > 0 ){
	key = ";";
	tmp1 = cut_str(tmp,key);
	chr2 = tmp1;
      }
      else{
	chr2 = chr1;
      }
    }

    if ( pos2 == -1 ){
      tmp = info ;
      tmp1 = cut_str(tmp,pos2_key);
      if ( tmp.size() > 0 ){
	pos2 = std::stol(tmp,&sz);
      }
    }

    if ( pos2 == -1 ){
      tmp = info ;
      tmp1 = cut_str(tmp,"END=");
      if ( tmp.size() > 0 ){
	pos2 = std::stol(tmp,&sz);
      }
    }

    /////////////////////////////////////////////////////////////////
    //
    // Reading STRAND
    //

    strand = "";
  
    key = "SVTYPE=BND";
    if ( info.find(key) != std::string::npos ){
      if ( alt[alt.size()-1] == '[' ) strand = "FR"; // fixed
      else if ( alt[alt.size()-1] == ']' ) strand = "FF"; // fixed
      else if ( alt[0] == '[' ) strand = "RR"; // fixed
      else if ( alt[0] == ']' ) strand ="RF"; // fixed
    }
  
    if ( alt[alt.size()-1] == '.' ){
      chr2="";
      pos2=-1;
      strand ="F";
    }
    else if ( alt[0] == '.' ){
      chr2="";
      pos2=-1;
      strand ="R";
    }

    if (strand.size() == 0){
      key="STRANDS=";
      tmp=info;
      cut_str ( tmp , key );
      if ( tmp.size() != 0 ){
	tmp = tmp.substr(0,2);
	if ( tmp == "+-" ) strand = "FR"; // fixed
	else if ( tmp == "++" ) strand = "FF"; // fixed
	else if ( tmp == "--" ) strand = "RR"; // fixed
	else if ( tmp == "-+" ) strand = "RF"; // fixed
      }
      else{
	key="CT=";
	tmp=info;
	cut_str(tmp,key);
	tmp=tmp.substr(0,4);
	if ( tmp == "5to3" ){
	  strand = "RF"; // fixed
	}
	else if ( tmp == "5to5" ){ 
	  strand = "RR"; // fixed
	}
	else if ( tmp == "3to3" ){
	  strand = "FF"; // fixed
	}
	else if ( tmp == "3to5" ){
	  strand = "FR"; // fixed
	}
      }
    }

    return line;
  }

  std::string to_string(){
    std::string tmp;
    tmp = chr1 + "\t" + std::to_string(pos1) + "\t" + sv_id + "\t" + ref + "\t" + alt + "\t" + qual + "\t" + filter + "\t" + info ;
    return tmp;
  }

};

using vcf_record_map=std::map < Position , std::vector < vcf_record > >;

class vcf_container{
private:

public:
  vcf_container(){
  }

  vcf_container(const std::string infile){
    read_vcf_file(infile);
  }

  ~vcf_container(){
  }

  vcf_record_map vcf_map;
  std::string vcf_file;
  std::string metainfo;
  std::string header;
  std::map < std::string , int > id_ref_map;
  std::map < int , std::string > ref_id_map;


  void read_vcf_file(const std::string infile){
    std::string tmp;
    vcf_record vcf;
    std::ifstream fin ( infile );
    while ( std::getline ( fin , tmp ) ){
      if ( tmp[1] == '#' ){
	metainfo += tmp + "\n";
      }
      else{
	header = tmp ;
	break;
      }
    }

    metainfo.pop_back();
  
    while ( std::getline ( fin , tmp ) ){
      vcf.parse(tmp);
      insert(vcf);
    }
  
    fin.close();
  }



  void insert(vcf_record vcf){
    Position Pos1;
    if ( id_ref_map.find(vcf.chr1) == id_ref_map.end()){
      int Size = id_ref_map.size();
      id_ref_map[vcf.chr1] = Size;
      ref_id_map[Size] = vcf.chr1;
    }
    Pos1.first = id_ref_map[vcf.chr1];
    Pos1.second= vcf.pos1;
    vcf_map[Pos1].push_back(vcf);
  }

  vcf_record_map::iterator begin(){return vcf_map.begin();};
  vcf_record_map::iterator end(){return vcf_map.end();};
  vcf_record_map::iterator find(Position Pos){
    return vcf_map.find(Pos);
  }


  std::vector < vcf_record > & operator [](Position Pos){
    return vcf_map[Pos];
  }






  void write(){
    std::cout << metainfo << "\n";
    std::cout << header << "\n";
    for ( auto & i : vcf_map ){
      for ( auto & j : i.second ){
	std::cout << j.to_string() << "\n";
      }
    }
  }



  void fwrite(std::string outfile){
    std::ofstream fout (outfile.c_str());
    fout << metainfo << "\n";
    fout << header << "\n";
    for ( auto & i : vcf_map ){
      for ( auto & j : i.second ){
	fout << j.to_string() << "\n";
      }
    }

    fout.close();
  }



  std::size_t size(){
    std::size_t count = 0 ;
    for ( auto & i : vcf_map ){
      count += i.second.size();
    }
    return count;
  }
 

};


void build_ref_id_map ( vcf_container & container, std::string infile){
  vcf_record vcf_line;
  std::string tmp;
  std::string id;
  int count = 0 ;

  std::ifstream fin;
  const std::string contigidkey="##contig=<ID=";
  fin.open( infile.c_str() );
  while ( std::getline ( fin , tmp ) ){
    if ( tmp.find(contigidkey) != std::string::npos ){
      tmp = tmp.substr(contigidkey.size());
      std::size_t found = tmp.find(",");
      id = tmp.substr(0,found);
      container.id_ref_map[id] = count;
      container.ref_id_map[count] = id;
      count ++;
    }
  }

  container.id_ref_map["."]=-1;
  container.ref_id_map[-1]=".";

}



/////////////////////////////////////////////////////////////////////////////////////
//
// Main
//


int main ( int argc , char ** argv ){
  if ( argc ==1 ){
    std::cout << "Usage:  remove_overlapping_sv  tumor.vcf  normal.vcf >  output.vcf\n";
    return 0;
  }

  double cputime0 = get_cpu_time();
  double walltime0 = get_wall_time();

  std::string tmp;
  vcf_record vcf;
  vcf_container container_t;
  vcf_container container_n;
  vcf_container container;

  std::string tumor_vcf = argv[1];
  std::string normal_vcf= argv[2];

  build_ref_id_map ( container_t , tumor_vcf );

  container_n=container_t;

  /////////////////////////////////////////////////////////////////////////////////////
  //
  // Reading tumor vcf
  //

  std::cerr << "Reading " << tumor_vcf << "\n";
  // container_t.read_vcf_file(tumor_vcf);
  std::ifstream fin ( tumor_vcf);
  while ( std::getline ( fin , tmp ) ){
    std::cout << tmp << "\n";
    if ( tmp[1] != '#' ){
      break;
    }
  }
  while ( std::getline ( fin , tmp ) ){
    vcf.parse(tmp);
    container_t.insert(vcf);
  }
  fin.close();

  /////////////////////////////////////////////////////////////////////////////////////
  //
  // Reading normal vcf
  //

  std::cerr << "Reading " << normal_vcf << "\n";
  // container_n.read_vcf_file(normal_vcf);
  fin.open ( normal_vcf );

  while ( std::getline ( fin , tmp ) ){
    if ( tmp[1] != '#' ){
      break;
    }
  }
  while ( std::getline ( fin , tmp ) ){
    vcf.parse(tmp);
    container_n.insert(vcf);
  }
  fin.close();

  
  std::map < int , std::map < int , int > > control_pos;

  for ( auto & i : container_n ){
    control_pos[i.first.first][i.first.second]=1;
    for ( auto & j : i.second ) {
      control_pos[container_n.id_ref_map[j.chr2]][j.pos2]=1;
    }
  }
  

  /////////////////////////////////////////////////////////////////////////////////////
  //
  // Remove germline
  //
  
  std::cerr << "Removing germline SVs\n";
  //container = remove_germline ( container_t, container_n );
  Position Pos1;
  Position Pos2;
  // std::string strand;

  Position Pos2_n;
  // std::string strand_n;

  for ( auto & i : container_t ){
    bool check = 1 ;
    for ( int x = -10 ; x <= 10 ; x ++ ){
      Pos1 = i.first;
      Pos1.second += x;

      if ( control_pos[Pos1.first].find(Pos1.second) != control_pos[Pos1.first].end() ){
	check = 0 ;
	// break;
      }
    }
    if ( check == 1 ){
      for ( auto & vcf : i.second ){
	std::cout << vcf.to_string() << "\n";
	// container.insert(vcf);
      }
    }
    else{
      for ( auto & vcf : i.second ){
	bool check = true ;

	for ( int x = -10 ; x <= 10 ; x ++ ){
	  Pos2.first = container_n.id_ref_map[vcf.chr2];
	  Pos2.second = vcf.pos2 + x;

	  if ( control_pos[Pos2.first].find(Pos2.second) != control_pos[Pos2.first].end() ){
	    check = 0;
	    // break;
	  }
	}

	if ( check == 1 ){
	  std::cout << vcf.to_string() << "\n";
	  // container.insert(vcf);
	}
      }
    }
  }


  

  //////////////////////////////////////////////////////////////////////////////////////////////
  //
  // End of program
  //

  double walltime1 = get_wall_time();
  double cputime1 = get_cpu_time();
  auto end = std::chrono::system_clock::now();
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);

  std::cerr << "===============================================================\n";
  std::cerr << "Total CPU time: " << cputime1 - cputime0 << " sec\n";
  std::cerr << "Total wall-clock time: " << walltime1 - walltime0 << " sec\n";
  std::cerr << "Peak memory: " << std::to_string((double)get_peak_memory()/(1024*1024)) << " GB\n";
  std::cerr << std::ctime(&end_time) << "\n";
  std::cerr << "===============================================================\n";



  return 0;
}



#endif
