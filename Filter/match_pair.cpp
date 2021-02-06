#include <string>
#include <fstream>
#include <iostream>
#include <unordered_set>

int main ( int argc , char ** argv){
  if ( argc == 0 ){
    std::cout << "Usage: match_pair input_1.fastq intput_2.fastq\n";
    return 0;
  }

  std::string id1;
  std::string id2;

  std::string seq1;
  std::string seq2;

  std::string desc1;
  std::string desc2;

  std::string qual1;
  std::string qual2;

  std::string infile1(argv[1]);
  std::string infile2(argv[2]);

  std::unordered_set < std::string > id1_set;
  std::unordered_set < std::string > id2_set;
  std::unordered_set < std::string > inter_id_set;


  std::ifstream fin1 ( infile1.c_str() );
  std::ifstream fin2 ( infile2.c_str() );

  while ( std::getline ( fin1 , id1 ) ){
    std::getline ( fin1 , seq1);
    std::getline ( fin1 , desc1);
    std::getline ( fin1 , qual1);

    id1_set.insert ( id1 );
  }

  while ( std::getline ( fin2 , id2 ) ){
    std::getline ( fin2 , seq2);
    std::getline ( fin2 , desc2);
    std::getline ( fin2 , qual2);

    id2_set.insert ( id2 );
  }

  fin1.close();
  fin2.close();

  for ( auto i : id1_set ){
    if ( id2_set.find(i) != id2_set.end() ){
      inter_id_set.insert(i);
    }
  }

  id1_set.clear();
  id2_set.clear();

  fin1.open ( infile1.c_str() );
  fin2.open ( infile2.c_str() );

  std::string outfile1 = infile1 + "_idmatched";
  std::string outfile2 = infile2 + "_idmatched";

  std::ofstream fout1 ( outfile1.c_str() );
  std::ofstream fout2 ( outfile2.c_str() );


  while ( std::getline ( fin1 , id1 ) ){
    std::getline ( fin1 , seq1);
    std::getline ( fin1 , desc1);
    std::getline ( fin1 , qual1);

    if ( inter_id_set.find(id1) != inter_id_set.end() ){
      fout1 << id1 << "\n" << seq1 << "\n" << desc1 << "\n" << qual1 << "\n";
    }
  }
  
  while ( std::getline ( fin2 , id2 ) ){
    std::getline ( fin2 , seq2);
    std::getline ( fin2 , desc2);
    std::getline ( fin2 , qual2);

    if ( inter_id_set.find(id2) != inter_id_set.end() ){
      fout2 << id2 << "\n" << seq2 << "\n" << desc2 << "\n" << qual2 << "\n";
    }
  }
  
  fin1.close();
  fin2.close();


  return 0;
}
