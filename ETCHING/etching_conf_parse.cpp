//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#include <string>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <stdio.h>
#include <limits.h>

void parse_usage(){
  std::cout << "Usage: etching_conf_parse input.conf\n";
}

int main ( int argc , char ** argv ){
  if ( argc == 1 ) {
    parse_usage();
    return 0;
  }
  
  std::string output;

  std::string tmp;
  std::string key;
  std::string value;
  std::ifstream fin ( argv [1] );
  std::string tmpfile = argv[1] ;
  std::string command;
  tmpfile += ".tmp";
  while ( std::getline ( fin , tmp ) ){
    if ( tmp[0] == '#' || tmp.size() == 0 ){
    }
    else{
      for ( unsigned int i = 0 ; i < tmp.size() ; i ++ ){
	if (isspace(tmp[i])){
	  key = tmp.substr(0,i);
	  tmp = tmp.substr(i);
	  for ( unsigned int j = 0 ; j < tmp.size() ; j ++ ){
	    if (!isspace(tmp[j])){
	      value = tmp.substr(j);
	      break;
	    }
	  }
	  break;
	}
      }

      // if ( key == "sample" || key == "filter_fq" || key == "filter_fa" || key == "annotation" ){
      if ( key == "sample" || key == "filter_fq" || key == "filter_fa" || key == "annotation" || key == "filter_db"){
	if ( value[0] != '/' && value.substr(0,2) != "~/" && value[0] != '$' ){
	  char cwd[PATH_MAX];
	  getcwd(cwd, sizeof(cwd));
	  std::string path=cwd;
	  value = path + "/" + value;
	  std::cout << key << "\t" << value << "\n";
	}
	else {
	  std::cout << key << "\t" << value << "\n";
	}
      }
      // else if ( key == "filter_db" ){
      // 	std::cout << key << "\t" << value << "\n";	
      // }
      else if ( key == "genome" ){
	if ( value[0] != '/' && value.substr(0,2) != "~/" && value[0] != '$' ){
	  char cwd[PATH_MAX];
	  getcwd(cwd, sizeof(cwd));
	  std::string path=cwd;
	  value = path + "/" + value;
	  std::cout << key << "\t" << value << "\n";
	}
	else{
	  std::cout << key << "\t" << value << "\n";
	}
      }
      else {
	std::cout << key << "\t" << value << "\n";
      }
      
    }
  }
  fin.close();

  return 0;
}
