#include <api/BamAlignment.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>

#include "Peak_memory_measure.hpp"
#include "CPU_TIME_measure.hpp"

#include "my_vcf.hpp"

int main ( int argc , char ** argv ){
  if ( argc !=4 ){
    std::cout << "Usage: extract_BP_read  input.vcf  input.bam  output.fastq\n";
    return -1;
  }

  std::string input_vcf = argv[1];
  std::string input_bam = argv[2];
  std::string output = argv[3];

  BamTools::BamReader reader;
  BamTools::RefVector references;
  BamTools::BamAlignment al;

  /////////////////////////////////////////////////////////////////////////////
  //
  // Bam initializing
  //
  
  reader.Open(input_bam);
  if ( !reader.OpenIndex(input_bam + ".bai") ) {
    reader.CreateIndex();
  }
  references = reader.GetReferenceData();


  /////////////////////////////////////////////////////////////////////////////
  //
  // Read vcf
  //

  std::map < std::string , int > id_ref_map;
  std::map < int , std::string > ref_id_map;

  for ( std::size_t i = 0 ; i < references.size() ; i ++ ){
    id_ref_map[references[i].RefName] = i ;
    ref_id_map[i] = references[i].RefName ;
  }
  id_ref_map["."]=-1;
  ref_id_map[-1]=".";

  std::map < int , std::map < int , int > > BP_position_map;

  std::ifstream fin ( input_vcf );
  std::string tmp;
  std::string chr_str;
  int chr;
  int pos;
  while ( std::getline ( fin , tmp ) ){
    if ( tmp[1] != '#' ){
      break;
    }
  }
  while ( std::getline ( fin , tmp ) ){
    std::stringstream ss ( tmp );
    ss >> chr_str >> pos ;
    chr = id_ref_map[chr_str];
    BP_position_map[chr][pos]=1;
  }

  fin.close();

  /////////////////////////////////////////////////////////////////////////////
  //
  // Scan bam
  //
  int count (0);
  int chr_num=-1;
  std::ofstream fout ( output.c_str() );
  while ( reader.GetNextAlignment(al) ){
    if (chr_num != al.RefID ){
      std::cerr << ref_id_map[al.RefID] << "\n";
      chr_num = al.RefID;
    }
    int chr1 = al.RefID;
    int pos1 = al.Position;
    int pos2 = al.GetEndPosition() - 1;
    
    for ( auto & i : BP_position_map[chr1] ){
      int pos = i.first;
      if ( pos1 - 100 < pos && pos < pos2 + 100 ){
	fout << "@" << count ++ << "\n"
	     << al.QueryBases << "\n"
	     << "+\n"
	     << al.Qualities << "\n";
	break;
      }
    }
    
  }
  fout.close();
  return 0;
}
