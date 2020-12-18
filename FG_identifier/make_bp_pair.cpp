//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#include "make_bp_pair.hpp"

int make_bp_pair(std::string input,std::string output){
  std::ifstream fin(input.c_str());
  std::ofstream fout(output.c_str());

  std::string chr1;
  std::string chr2;
  int pos1;
  int pos2;
  std::string id;
  std::string ref;
  std::string alt;
  std::string qual;
  std::string filter;
  std::string info;

  std::string tmp;

  const std::string START="#CHROM";
  std::string in_line;
  
  std::size_t found;

  while ( getline ( fin , in_line ) ){
    if ( in_line.find(START) != std::string::npos ){
      break;
    }
  }


  while ( fin >> chr1 >> pos1 >> id >> ref >> alt >> qual >> filter ){
    getline ( fin , info );
    // FF
    if ( alt.back() == '[' ){
      found = alt.find("[");
      tmp = alt.substr(found+1);
      found = tmp.find("[");
      tmp = tmp.substr(0,found);
      found = tmp.find(":");
      chr2=tmp.substr(0,found);
      pos2=atoi(tmp.substr(found+1).c_str());
      fout << chr1 << "\t" << pos1 << "\t" << chr2 << "\t" << pos2 << "\tFF\n";
    }
      
    // FR
    else if ( alt.back() == ']' ){
      found = alt.find("]");
      tmp = alt.substr(found+1);
      found = tmp.find("]");
      tmp = tmp.substr(0,found);
      found = tmp.find(":");
      chr2=tmp.substr(0,found);
      pos2=atoi(tmp.substr(found+1).c_str());
      fout << chr1 << "\t" << pos1 << "\t" << chr2 << "\t" << pos2 << "\tFR\n";
    }
      
    // RF
    else if ( alt.front() == '[' ){
      found = alt.find("[");
      tmp = alt.substr(found+1);
      found = tmp.find("[");
      tmp = tmp.substr(0,found);
      found = tmp.find(":");
      chr2=tmp.substr(0,found);
      pos2=atoi(tmp.substr(found+1).c_str());
      fout << chr1 << "\t" << pos1 << "\t" << chr2 << "\t" << pos2 << "\tRF\n";
    }
      
    // RR
    else if ( alt.front() == ']' ){
      found = alt.find("]");
      tmp = alt.substr(found+1);
      found = tmp.find("]");
      tmp = tmp.substr(0,found);
      found = tmp.find(":");
      chr2=tmp.substr(0,found);
      pos2=atoi(tmp.substr(found+1).c_str());
      fout << chr1 << "\t" << pos1 << "\t" << chr2 << "\t" << pos2 << "\tRR\n";
    }
  
    else if ( info.find("SVTYPE=DEL") != std::string::npos ){
      found = info.find(";END=");
      if ( found == std::string::npos ) found = info.find("\tEND=");
      tmp = info.substr(found+5);
      found = tmp.find(";");
      tmp = tmp.substr(0,found);
      pos2=atoi(tmp.substr(0,found).c_str());
      chr2=chr1;
      fout << chr1 << "\t" << pos1 << "\t" << chr2 << "\t" << pos2 << "\tFF\n";
      fout << chr2 << "\t" << pos2 << "\t" << chr1 << "\t" << pos1 << "\tRR\n";
    }
    else if ( info.find("SVTYPE=DUP") != std::string::npos ){
      found = info.find(";END=");
      if ( found == std::string::npos ) found = info.find("\tEND=");
      tmp = info.substr(found+5);
      found = tmp.find(";");
      tmp = tmp.substr(0,found);
      pos2=atoi(tmp.substr(0,found).c_str());
      chr2=chr1;
      fout << chr1 << "\t" << pos1 << "\t" << chr2 << "\t" << pos2 << "\tRR\n";
      fout << chr2 << "\t" << pos2 << "\t" << chr1 << "\t" << pos1 << "\tFF\n";
    }
    else if ( info.find("SVTYPE=INV") != std::string::npos ){
      found = info.find(";END=");
      if ( found == std::string::npos ) found = info.find("\tEND=");
      tmp = info.substr(found+5);
      found = tmp.find(";");
      tmp = tmp.substr(0,found);
      pos2=atoi(tmp.substr(0,found).c_str());
      chr2=chr1;
      fout << chr1 << "\t" << pos1 << "\t" << chr2 << "\t" << pos2 << "\tFR\n";
      fout << chr1 << "\t" << pos1 << "\t" << chr2 << "\t" << pos2 << "\tRF\n";
      fout << chr2 << "\t" << pos2 << "\t" << chr1 << "\t" << pos1 << "\tFR\n";
      fout << chr2 << "\t" << pos2 << "\t" << chr1 << "\t" << pos1 << "\tRF\n";
    }
    else if ( info.find("SVTYPE=TRA") != std::string::npos ){
      found = info.find(";END=");
      if ( found == std::string::npos ) found = info.find("\tEND=");
      tmp = info.substr(found+5);
      found = tmp.find(";");
      tmp = tmp.substr(0,found);
      pos2=atoi(tmp.substr(0,found).c_str());
      chr2=chr1;
      if ( info.find("CT=5to3") != std::string::npos ){ // FF
	fout << chr1 << "\t" << pos1 << "\t" << chr2 << "\t" << pos2 << "\tFF\n";
	fout << chr2 << "\t" << pos2 << "\t" << chr1 << "\t" << pos1 << "\tRR\n";      
      } 
      else if ( info.find("CT=3to3") != std::string::npos ){ //RF
	fout << chr1 << "\t" << pos1 << "\t" << chr2 << "\t" << pos2 << "\tRF\n";
	fout << chr2 << "\t" << pos2 << "\t" << chr1 << "\t" << pos1 << "\tRF\n";      
      }
      else if ( info.find("CT=5to5") != std::string::npos ){ // FR
	fout << chr1 << "\t" << pos1 << "\t" << chr2 << "\t" << pos2 << "\tFR\n";
	fout << chr2 << "\t" << pos2 << "\t" << chr1 << "\t" << pos1 << "\tFR\n";      
      }
      else if ( info.find("CT=3to5") != std::string::npos ){ // RR
	fout << chr1 << "\t" << pos1 << "\t" << chr2 << "\t" << pos2 << "\tRR\n";
	fout << chr2 << "\t" << pos2 << "\t" << chr1 << "\t" << pos1 << "\tFF\n";      
      }
    }
  }

  fout.close();
  fin.close();

  return 0;
}
