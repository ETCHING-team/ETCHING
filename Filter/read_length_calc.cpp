#include <string>
#include <iostream>
#include <gzstream.h>
#include <functional>
#include <vector>

#include <api/BamAlignment.h>
#include <api/BamReader.h>


std::vector <std::string> fastq_count ( std::string infile , int max_count){
  std::vector < std::string > output; 
  std::ifstream fin (infile.c_str());
  std::string id;
  std::string seq;
  std::string desc;
  std::string qual;

  int count=0 ;

  while ( getline ( fin , id ) && getline ( fin , seq ) && getline ( fin , desc ) && getline ( fin , qual ) && count < max_count){
    output.push_back ( seq );
    count ++;
  }
  
  fin.close();
  return output;
}

std::vector <std::string> fastq_gz_count ( std::string infile , int max_count){
  std::vector < std::string > output;
  igzstream fin;
  fin.open(infile.c_str());
  std::string id;
  std::string seq;
  std::string desc;
  std::string qual;

  int count=0 ;

  while ( getline ( fin , id ) && getline ( fin , seq ) && getline ( fin , desc ) && getline ( fin , qual ) && count < max_count){
    output.push_back ( seq );
    count ++;
  }
  
  fin.close();
  return output;
}


std::vector <std::string> bam_count ( std::string infile , int max_count){
  std::vector < std::string > output;

  BamTools::BamReader reader;
  BamTools::BamAlignment al;

  reader.Open ( infile );
  
  int count = 0 ;

  while ( reader.GetNextAlignment(al) && count < max_count ){
    output.push_back(al.QueryBases);
    count ++;
  }

  return output;

}


int main ( int argc , char ** argv ){
  if ( argc != 3 ){
    std::cout << "Usage: read_length_calc input max_count\n";
    return 0;
  }

  std::string infile = argv[1];
  int max_count = atoi ( argv[2] );

  int read_length_sum =0 ;
  int read_length_average =0 ;
  int count = 0;

  std::string id;
  std::string seq;
  std::string desc;
  std::string qual;

  std::vector < std::string > (*fcnPtr)(std::string, int);
  std::vector < std::string > seq_vec;
  
  int Size = infile.size();
  if ( infile.substr(Size-6) == ".fastq" ||
       infile.substr(Size-3) == ".fq" ){
    fcnPtr = fastq_count;
  }
  else if ( infile.substr(Size-9) == ".fastq.gz" ||
	    infile.substr(Size-6) == ".fq.gz" ){
    fcnPtr = fastq_gz_count;
  }
  else if ( infile.substr(Size-4) == ".bam" ){
    fcnPtr = bam_count;
  }

  seq_vec = fcnPtr ( infile , max_count );
  count = seq_vec.size();

  for ( auto i : seq_vec ) {
    read_length_sum += i.size();
    //std::cout << i.size() << "\n";
  }

  read_length_average = read_length_sum / count ;
  std::cout << read_length_average << "\n";

  return 0;
}
