#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

int main ( int argc , char ** argv ){
  if ( argc == 1 ){
    std::cout << "Usage: cut_by_score etching.vcf cutoff(0-1)\n";
    return 0;
  }

  std::string infile ( argv[1] );
  std::ifstream fin ( infile.c_str() );

  std::size_t sz;
  double cutoff = std::stod ( argv[2] , &sz );

  std::string line;
  while ( std::getline ( fin , line ) ){
    std::cout << line << "\n";
    if ( line[0] == '#' && line[1] != '#' ){
      break;
    }
  }

  const std::string score_key="SVSCORE";

  while ( std::getline ( fin , line ) ){
    
    std::stringstream line_ss ( line );
    std::string info;

    double score;

    for ( std::size_t i = 0 ; i < 8 ; i ++ ){
      if ( i == 5 ){
	line_ss >> score;
      }
      else{
	line_ss >> info ;
      }
    }
    
    if ( line.find(score_key) != std::string::npos ){
      std::stringstream info_ss ( info );
      std::string info_elem;
      while ( std::getline ( info_ss , info_elem , ';') ){

	std::stringstream key_ss ( info_elem ); 
	std::string key;
	std::string value;
      
	std::getline ( key_ss , key , '=' );
	
	if ( key == score_key ){
	  std::getline ( key_ss , value , '=' );
	  score = std::stod ( value , &sz );
	  break;
	}
      }
    }
    
    if ( score >= cutoff ){
      std::string lowqual = "LOWQUAL";
      std::size_t found = line.find ( lowqual );
      if ( found != std::string::npos ){
	std::string pass = "PASS";
	line.replace ( found, lowqual.size(), pass );
      }
      std::cout << line << "\n";
    }
  }

  fin.close();
  return 0;
}
