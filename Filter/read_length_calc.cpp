#include <string>
#include <iostream>
#include <gzstream.h>
#include <functional>
#include <vector>

#include <api/BamAlignment.h>
#include <api/BamReader.h>


std::vector <std::string> fastq_count ( std::string infile , int64_t max_count){
  std::vector < std::string > output; 
  std::ifstream fin (infile.c_str());
  std::string id;
  std::string seq;
  std::string desc;
  std::string qual;

  int64_t count=0 ;

  while ( getline ( fin , id ) && getline ( fin , seq ) && getline ( fin , desc ) && getline ( fin , qual ) && count < max_count){
    output.push_back ( seq );
    count ++;
  }
  
  fin.close();
  return output;
}

std::vector <std::string> fastq_gz_count ( std::string infile , int64_t max_count){
  std::vector < std::string > output;
  igzstream fin;
  fin.open(infile.c_str());
  std::string id;
  std::string seq;
  std::string desc;
  std::string qual;

  int64_t count=0 ;

  while ( getline ( fin , id ) && getline ( fin , seq ) && getline ( fin , desc ) && getline ( fin , qual ) && count < max_count){
    output.push_back ( seq );
    count ++;
  }
  
  fin.close();
  return output;
}


std::vector <std::string> bam_count ( std::string infile , int64_t max_count){
  std::vector < std::string > output;

  BamTools::BamReader reader;
  BamTools::BamAlignment al;

  reader.Open ( infile );
  
  int64_t count = 0 ;

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
  int64_t max_count = atoi ( argv[2] );

  int64_t read_length_sum =0 ;
  int64_t read_length_average =0 ;
  int64_t count = 0;

  std::string id;
  std::string seq;
  std::string desc;
  std::string qual;

  std::vector < std::string > (*fcnPtr)( std::string, int64_t );
  std::vector < std::string > seq_vec;

  int64_t Size = infile.size();

  int64_t size9 = Size - 9;
  if (size9 < 0 ) size9=0;

  int64_t size6 = Size - 6;
  if (size6 < 0 ) size6=0;

  int64_t size4 = Size - 4;
  if (size4 < 0 ) size4=0;

  int64_t size3 = Size - 3;
  if (size3 < 0 ) size3=0;

  if ( infile.substr(size6) == ".fastq" ||
       infile.substr(size3) == ".fq" ){
    fcnPtr = fastq_count;
  }
  else if ( infile.substr(size9) == ".fastq.gz" ||
	    infile.substr(size6) == ".fq.gz" ){
    fcnPtr = fastq_gz_count;
  }
  else if ( infile.substr(size4) == ".bam" ){
    fcnPtr = bam_count;
  }

  seq_vec = fcnPtr ( infile , max_count );
  count = seq_vec.size();

  for ( auto i : seq_vec ) {
    read_length_sum += i.size();
  }

  read_length_average = read_length_sum / count ;
  std::cout << read_length_average << "\n";

  return 0;
}
