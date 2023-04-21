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
#include <sstream>

void typer_usage(){
  std::cout << "Usage:  etching_typer  input.vcf  genome.fa.fai\n" ;
}


int main ( int argc , char ** argv ){
  if ( argc != 3 ){
  typer_usage();
  return 0;
  }
  
  std::string infile = argv[1];
  std::string index = argv[2];

  std::string id;
  int length;

  std::string tmp;
  std::vector < std::string > metainfo;
  VCF_LINE vcf_line;
  VCF container;

  int count = 0 ;
  
  std::ifstream fin;
  fin.open( index.c_str() );
  while ( fin >> id ){
    fin >> length;
    container.id_ref_map[id] = count;
    container.ref_id_map[count] = id;
    count ++;
  }

  fin.close();
  fin.open( infile.c_str() );

  /////////////////////////////////////////////////////////////////
  //
  // Reading metainfo
  //

  //bool is_etching = 0 ;

  while ( std::getline ( fin , tmp ) ){
    if ( tmp.substr(0,2) == "##" ){
      metainfo.push_back(tmp);
      //if ( tmp.find("##source=ETCHING") != 0 ) {
      //is_etching = 1;
      //}
    }
    else break;
  }
  metainfo.push_back(tmp);

  
  /////////////////////////////////////////////////////////////////
  //
  // fill VCF container
  //

  while ( std::getline ( fin , tmp ) ){
    //if ( is_etching) vcf_line.parse_etching(tmp);
    //else 
    vcf_line.parse(tmp);
    container.insert(vcf_line);
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
  // Typing
  //
  std::cout << container.metainfo << "\n";
  std::cout << container.header << "\n";
  for ( auto & i : container.vcf_map ){
    std::string svtype;
    if ( i.second.chr1 == i.second.chr2 ){
      if ( i.second.pos1 < i.second.pos2 ){
	if ( i.second.strand == "FR" ) svtype = "DEL";
	else if ( i.second.strand == "RF" ) svtype = "DUP";
	else if ( i.second.strand == "FF" ) svtype = "INV";
	else if ( i.second.strand == "RR" ) svtype = "INV";
      }
      if ( svtype.size() != 0 ){
	std::stringstream additional_info(i.second.info);
	std::vector < std::string > additional_info_vec;
	std::string token;
	while ( additional_info >> token ) additional_info_vec.push_back(token);
	std::string INFO = additional_info_vec[0];
	std::stringstream info_stream(INFO);
	std::vector < std::string > INFO_vec;
	std::string info;
	while ( std::getline ( info_stream , info , ';' ) ){
	  if ( info.substr(0,7) == "SVTYPE=" ){
	    info="SVTYPE="+svtype;
	    i.second.alt="<"+svtype+">";
	  }
	  //if ( info.substr(0,7) == "MATEID=" ){
	  //info="SOURCEBND=" + i.second.sv_id + "," + info.substr(7) ;
	  //}
	  INFO_vec.push_back(info);
	}
	additional_info_vec[0].clear();
	for ( auto j : INFO_vec ){
	  additional_info_vec[0] += j + ";";
	}
	additional_info_vec[0].pop_back();
	i.second.info.clear();
	for ( auto j : additional_info_vec ){
	  i.second.info += j + "\t";
	}
	i.second.info.pop_back();
	std::cout << i.second.to_string() << "\n";
      }
    }
    else{
      std::cout << i.second.to_string() << "\n";
    }
  }
  
  return 0;
}
