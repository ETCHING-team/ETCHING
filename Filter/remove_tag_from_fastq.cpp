#include <string>
#include <iostream>
#include <fstream>

bool tag_check ( std::string id ){
  std::string pref = id.substr( id.size() - 2);
  return ( pref == "/1" || pref == "/2" ? 1 : 0 );
}

int main ( int argc , char ** argv ){
  if ( argc == 1 ){
    std::cout << "remove_tag_from_fastq fastq\n";
    return 0;
  }

  std::string infile = argv[1];

  std::string id;
  std::string seq;
  std::string desc;
  std::string qual;

  std::ifstream fin ( infile.c_str() );
  while ( fin >> id >> seq >> desc >> qual ){
    if ( tag_check ( id ) ){
      id = id.substr ( 0, id.size() - 2 );
    }
    std::cout << id << "\n" << seq << "\n" << desc << "\n" << qual << "\n";
  }
  fin.close();

  return 0;
}
