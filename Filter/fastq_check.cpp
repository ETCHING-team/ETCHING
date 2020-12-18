#include <iostream>
#include <fstream>
#include <gzstream.h>
#include <string>

void usage(){
  std::cout << "fastq_check input.fastq(.gz)\n";
  std::cout << "\t" << "Return value: 0: No issues\n";
  std::cout << "\t" << "            : 1: tab in id\n";
  std::cout << "\t" << "            : 2: Broken (or not) fastq\n";

}

int main ( int argc , char ** argv ){
  if ( argc !=2 ){
    usage();
    return 0;
  }
  
  std::string infile = argv[1];

  std::string id;
  std::string seq;
  std::string desc;
  std::string qual;

  if ( infile.substr(infile.size() - 3) != ".gz" ){
    std::ifstream fin ( infile.c_str() );
    std::getline ( fin , id  );
    std::getline ( fin , seq );
    std::getline ( fin , desc);
    std::getline ( fin , qual);
    fin.close();
  } else {
    igzstream fin ;
    fin.open( infile.c_str() );
    std::getline ( fin , id  );
    std::getline ( fin , seq );
    std::getline ( fin , desc);
    std::getline ( fin , qual);
    //fin.Close();
  }
  
  for ( std::size_t i = 0 ; i < id.size() ; i ++ ){
    if ( id[i] == 9 ) {
      std::cout << "Warning!!! Input fastq file as \\t character(s).\n";
      return 1;
    }
  }

  if ( seq.size() != qual.size() || desc[0] != '+' || id[0] != '@' ){
    std::cout << "Warning!!! Input file is broken (or not) fastq.\n";
    return 2;
  }

  return 0;
}
