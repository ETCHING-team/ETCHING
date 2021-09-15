//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#include "etching_caller.hpp"
#include <unordered_map>
void initial_bp_call (  const std::string input_bam, const std::string prefix)
{

  BamTools::BamAlignment al;
  BamTools::BamReader reader;
  BamTools::SamHeader header;
  BamTools::RefVector references;

  reader.Open(input_bam);
  if ( !reader.OpenIndex(input_bam + ".bai") ) {
    std::cerr << "!!!ERROR!!! Indexed bam file is required!!!\t";
    exit (EXIT_FAILURE);
  }
  header = reader.GetHeader();
  references = reader.GetReferenceData();

  std::cerr << "[Step 1: Finding initial BPs by scanning bam file]\n";
  std::cerr << "Scanning " << input_bam << "\n";

  double cputime0 = get_cpu_time();
  double walltime0 = get_wall_time();

  std::unordered_map < std::string , int > mapping_count;
  std::unordered_set < std::string > selected_id_set;
  std::unordered_set < std::string > split_read_set;

  
  //////////////////////////////////////////////////////////////////////
  //
  // Count softclip
  //
  
  std::cerr << "Counting soft- and hard-clip\n" ;

  std::map < Position , int > softclip_map_for;
  std::map < Position , int > softclip_map_rev;
  Position Pos;
  while ( reader.GetNextAlignment(al) ) {
    if ( al.CigarData[0].Type == 'S' || al.CigarData[0].Type == 'H' ){
      Pos.first = al.RefID;
      Pos.second = al.Position;
      if ( Pos.first >= 0 && Pos.second >= 0 ){
	softclip_map_rev [ Pos ] --;
	selected_id_set.insert(al.Name);
      }
    }
    if ( al.CigarData[al.CigarData.size()-1].Type == 'S' || al.CigarData[al.CigarData.size()-1].Type == 'H' ){
      Pos.first = al.RefID;
      Pos.second = al.GetEndPosition() - 1;
      if ( Pos.first >= 0 && Pos.second >= 0 ){
	softclip_map_for [ Pos ] ++ ;
	selected_id_set.insert(al.Name);
      }
    }
  }
  
  //////////////////////////////////////////////////////////////////////
  //
  // Print simply predicted BPs using SD_cutoff
  //
  
  std::cerr << "Printing initial BPs\n" ;
  

  std::string outfile = prefix + ".initial_bp.txt";
  std::ofstream fout ( outfile.c_str() );

    
  for ( auto & i : softclip_map_for ){
    fout << references[i.first.first].RefName << "\t" << std::to_string(i.first.second+1) << "\t" << i.second << "\n";
  }

  for ( auto & i : softclip_map_rev ){
    fout << references[i.first.first].RefName << "\t" << std::to_string(i.first.second+1) << "\t" << i.second << "\n";
  }



  double cputime1 = get_cpu_time();
  double walltime1 = get_wall_time();

  std::cerr << "End of Step 1\n";
  std::cerr << "CPU time: " << cputime1 - cputime0 << " sec\n";
  std::cerr << "Wall-clock time: " << walltime1 - walltime0 << " sec\n\n";

  fout.close();
}
