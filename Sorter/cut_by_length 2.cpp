#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

int main ( int argc , char ** argv ){
  if ( argc == 1 ){
    std::cout << "Usage: cut_by_length input.vcf length_cutoff\n";
    return 0;
  }

  std::string infile ( argv[1] );
  std::ifstream fin ( infile.c_str() );

  int cutoff = atoi ( argv[2] );

  std::string line;
  while ( std::getline ( fin , line ) ){
    std::cout << line << "\n";
    if ( line[0] == '#' && line[1] != '#' ){
      break;
    }
  }

  while ( std::getline ( fin , line ) ){
    
    std::string chr1;
    std::string chr2;

    int pos1;
    int pos2;

    const std::string chr2_key="CHR2";
    const std::string pos2_key="END";

    std::stringstream line_ss ( line );
    std::string info;

    for ( std::size_t i = 0 ; i < 8 ; i ++ ){
      if ( i == 0 ){
	line_ss >> chr1;
      }
      else if ( i == 1 ){
	line_ss >> pos1;
      }
      else{
	line_ss >> info ;
      }
    }
    

    std::stringstream info_ss ( info );
    std::string info_elem;
    while ( std::getline ( info_ss , info_elem , ';') ){

      std::stringstream key_ss ( info_elem ); 
      std::string key;
      std::string value;
      
      std::getline ( key_ss , key , '=' );
	
      if ( key == chr2_key ){
	key_ss >> chr2;
      }
      else if ( key == pos2_key ){
	key_ss >> pos2;
      }
    }

    if ( chr1 == chr2 ){
      int diff = pos2 - pos1;
      if ( diff < 0 ) diff *= -1;

      if ( diff >= cutoff ){
	std::cout << line << "\n";
      }
    }
    else{
      std::cout << line << "\n";
    }
  }

  fin.close();
  return 0;
}
