//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#include <string>
#include <map>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

int main ( int argc , char ** argv ){
  if ( argc == 1 ){
    std::cout << "Usage:  sort_gtf  gtf  > output\n";
    return 0;
  }

  std::string comment;
  std::map < std::pair < std::size_t, std::size_t > , std::vector < std::string > > gtf_map;
  std::map < std::string , std::size_t > chr_map;

  std::string chr;
  std::size_t first;
  std::string second;
  std::string third;
  std::size_t start;
  std::string last;

  std::string input;

  std::size_t count = 0 ;

  std::ifstream fin ( argv[1] ) ;
  for ( int i = 0 ; i < 5 ; i ++ ){
    getline ( fin , last );
    comment += chr + last + '\n';
  }
  while ( fin >> chr ){
    getline ( fin , last );
    if ( chr_map.find(chr) == chr_map.end() ){
      chr_map[chr]=count++;
    }
  }
  fin.close();

  fin.open(argv[1]);
  for ( int i = 0 ; i < 5 ; i ++ ){
    getline ( fin , last );
  }
  while ( fin >> chr >> second >> third >> start ){
    getline ( fin , last );
    first = chr_map[chr];
    input = chr + "\t" 
      + second + "\t" 
      + third + "\t"
      + std::to_string(start) 
      + last;
    gtf_map[std::make_pair(first,start)].push_back(input);
  }
  fin.close();

  std::cout << comment ;
  for ( auto & i : gtf_map ){
    for ( auto & j : i.second ){
      std::cout << j << "\n";
    }
  }
  return 0;
}
