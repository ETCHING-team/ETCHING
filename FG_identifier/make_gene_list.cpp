//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------


#include "make_gene_list.hpp"


void parsing_gff_info(std::string input, std::map < std::string , std::string > & info ){
  std::size_t found = input.find(";");
  std::string output = input.substr(found+1);
  if ( found != std::string::npos ){
    std::string tmp = input.substr(0,found);
    found = tmp.find("=");
    info[tmp.substr(0,found)]=tmp.substr(found+1);
    parsing_gff_info(output,info);
  }
}




int make_gene_list ( std::string arg1, std::string arg2, std::string arg3, std::string arg4){
  
  std::string id;
  std::string start;
  std::string end;
  std::string strand;
  std::string type;
  std::string name;
  std::string tmp;
  
  const std::string infile=arg1;
  const std::string feature=arg2;
  const std::string key=arg3;
  std::ifstream fin ( arg1.c_str() );
  std::ofstream fout ( arg4.c_str() );

  std::size_t found;

  if ( infile.substr(infile.size()-4) == ".gtf" ){
    while ( fin >> id >> tmp >> type >> start >> end >> tmp >> strand) {
      getline ( fin , tmp );
      if ( type == feature ){
	found = tmp.find(key);
	tmp=tmp.substr(found+key.size()+2);
	found = tmp.find("\"");
	tmp=tmp.substr(0,found);
	fout << id << "\t" << start << "\t" << end << "\t" << strand << "\t" << tmp << "\n";
      }
    }
  }
  else if ( infile.substr(infile.size()-4) == ".gff" ){
    while ( fin >> id >> tmp >> type >> start >> end >> tmp >> strand >> tmp) {
      getline ( fin , tmp );
      if ( type == feature){
	std::map < std::string , std::string > info;
	if ( isspace(tmp[0]) ) tmp = tmp.substr(1);
	parsing_gff_info(tmp,info);
	fout << id << "\t" << start << "\t" << end << "\t" << strand << "\t" << info[key] << "\n";
      }
    }
  }
  else {
    std::cout << "ERROR!!! Input annotation file must be gtf or gff file\n";
  }

  return 0;
}

