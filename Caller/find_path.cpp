//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#include "etching_caller.hpp"
#include "functions.hpp"


void GetReferenceData(const std::string input_bam, 
		      BamTools::RefVector & references, 
		      std::unordered_map < std::string , int > & id_ref_map)
{
  BamTools::BamAlignment al;
  BamTools::BamReader reader;
  std::string chr;

  reader.Open(input_bam);
  if ( !reader.OpenIndex(input_bam + ".bai") ) {
    reader.CreateIndex();
    //std::cout << "!!!ERROR!!! BAM file was not indexed!!!\n";
    //exit (EXIT_FAILURE);
  }
  
  references = reader.GetReferenceData();

  for ( std::size_t i = 0 ; i < references.size() ; i ++ ){
    id_ref_map[references[i].RefName] = i ;
  }
}






void Find_BP_paths(std::string input_bam, 
		   std::unordered_map < std::string , int > id_ref_map, 
		   std::map < Position, std::map < Position , std::map < std::string , int > > > & BP_graph, 
		   const int insert_size, 
		   const bool scanall )
{

  std::cout << "Building BP-graph using split-reads\n";
    
  Position Pos;

  BamTools::BamReader Reader;
  Reader.Open(input_bam);
  Reader.Jump(0,0);

  std::string Tag;
    
  std::map < int , std::string > ref_id_map;
  for ( auto & i : id_ref_map ){
    ref_id_map[i.second]=i.first;
  }
  
  // int strand;
    
  Position Pos_read;  // read position
  Position Pos_mate;  // mate position
  Position Pos_split; // split position
  Position Pos_unmap;
  Pos_unmap.first = -1;
  Pos_unmap.second= -1;
    
  BamTools::BamAlignment Al;

  int exID = - 1 ;
  while ( Reader.GetNextAlignment(Al) ) {
    if ( Al.RefID != exID ) {
      std::cout << ref_id_map[Al.RefID] << "\n";
      exID = Al.RefID;
    }
    if ( Al.GetTag<std::string>("SA",Tag) ){
      std::vector < std::string > tag_vec;
      while ( Tag.size() > 0 ){
	std::string tmp = parse_SA_tag(Tag);
	tag_vec.push_back( tmp );
      }
      if ( Al.IsMateMapped () ) {
	if ( ( Al.CigarData[Al.CigarData.size()-1].Type == 'S' || Al.CigarData[Al.CigarData.size()-1].Type == 'H' ) 
	     && ( Al.CigarData[0].Type != 'S' && Al.CigarData[0].Type != 'H' ) ){
	  std::string strand = "F";
	  
	  Pos_read.first = Al.RefID ;
	  Pos_read.second = Al.GetEndPosition()-1;
	  for ( auto tag : tag_vec){
	    Pos_split = find_split(tag,id_ref_map);
	    Pos_mate.first = Al.MateRefID;
	    Pos_mate.second = Al.MatePosition;
	    if ( Pos_split.second != -1 ){
	      strand += mate_strand (tag);
	      if ( scanall ){
		BP_graph[Pos_read][Pos_split][strand]++; 
	      }
	      else if ( FindConcordant ( Pos_read , Pos_split , Pos_mate , insert_size ) ){
		BP_graph[Pos_read][Pos_split][strand]++;
	      }
	    }
	  }
	}
	else if ( ( Al.CigarData[0].Type == 'S' || Al.CigarData[0].Type == 'H' )
		  && ( Al.CigarData[Al.CigarData.size()-1].Type != 'S' && Al.CigarData[Al.CigarData.size()-1].Type != 'H' ) ) {
	  std::string strand = "R";
	  
	  Pos_read.first = Al.RefID ;
	  Pos_read.second = Al.Position;
	  for ( auto tag : tag_vec){
	    Pos_split = find_split(tag,id_ref_map);
	    Pos_mate.first = Al.MateRefID;
	    Pos_mate.second = Al.MatePosition;
	    if ( Pos_split.second != -1 ){
	      strand += mate_strand (tag);
	      if ( scanall ){
		BP_graph[Pos_read][Pos_split][strand]++;
	      }
	      else if ( FindConcordant ( Pos_read , Pos_split , Pos_mate , insert_size ) ){
		BP_graph[Pos_read][Pos_split][strand]++;
	      }
	    }
	  }
	}
      }
      // Scan single breakend
      else {
	if ( ( Al.CigarData[Al.CigarData.size()-1].Type == 'S' || Al.CigarData[Al.CigarData.size()-1].Type == 'H' ) 
	     && ( Al.CigarData[0].Type != 'S' && Al.CigarData[0].Type != 'H' ) ){
	  std::string strand = "F";
	  Pos.first = Al.RefID ;
	  Pos.second = Al.GetEndPosition()-1;
	  BP_graph[Pos][Pos_unmap][strand]++; 
	}
	else if ( ( Al.CigarData[0].Type == 'S' || Al.CigarData[0].Type == 'H' )
		  && ( Al.CigarData[Al.CigarData.size()-1].Type != 'S' && Al.CigarData[Al.CigarData.size()-1].Type != 'H' ) ){
	  std::string strand = "R";
	  Pos.first = Al.RefID ;
	  Pos.second = Al.Position;
	  BP_graph[Pos][Pos_unmap][strand]++; 
	}
      }
    }
  }
  Reader.Close();
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
void remove_germline(std::string normal_bam, 
		   std::unordered_map < std::string , int > id_ref_map,
		   std::map < Position, std::map < Position , std::map < std::string , int > > > & BP_graph,
		   const int insert_size, const bool scanall )
{
  std::cout << "Building BP-graph using split-reads\n";

  std::map < Position, std::map < Position , std::map < std::string , int > > > normal_BP_graph;
  std::map < Position, std::map < Position , std::map < std::string , int > > > output_BP_graph;
    
  Find_BP_paths(normal_bam, id_ref_map, normal_BP_graph, insert_size, scanall);

  Position Pos1;
  Position Pos2;
  std::string strand;
  int value;
  
  for ( auto & i : BP_graph ){
    Pos1 = i.first;
    for ( auto & j : i.second ){
      Pos2 = j.first;
      for ( auto & k : j.second ){
	strand = k.first;
	value = k.second;
	  
	if ( normal_BP_graph.find(Pos1) == normal_BP_graph.end() ){
	  output_BP_graph[Pos1][Pos2][strand] = value;
	}
	else if ( normal_BP_graph[Pos1].find(Pos2) == normal_BP_graph[Pos1].end() ) {
	  output_BP_graph[Pos1][Pos2][strand] = value;
	}
	else if ( normal_BP_graph[Pos1][Pos2].find(strand) == normal_BP_graph[Pos1][Pos2].end() ){
	  output_BP_graph[Pos1][Pos2][strand] = value;
	}
      }
    }
  }
  
  BP_graph.clear();
  BP_graph = output_BP_graph;
}
*/

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Writing_files( const std::string prefix, 
		    std::map < Position , std::map < Position , std::map < std::string , int > > > & BP_graph,
		    BamTools::RefVector references )
{
  std::string outfile_pair = prefix + ".pair.txt";
  std::ofstream fout;
  fout.open(outfile_pair.c_str());
  for ( auto & i : BP_graph ){
    for ( auto & j : i.second ){
      for ( auto & k : j.second ){
	if ( j.first.second > -1 ){
	  fout << references[i.first.first].RefName << "\t" << i.first.second + 1 << "\t"
	       << references[j.first.first].RefName << "\t" << j.first.second + 1 << "\t"
	       << k.first << "\t" << k.second << "\n";
	}
	else{
	  fout << references[i.first.first].RefName << "\t" << i.first.second + 1 << "\t"
	       << "." << "\t" << -1 << "\t"
	       << k.first << "\t" << k.second << "\n";
	}
      }
    }
  }
  fout.close();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  MAIN
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////



void find_path ( const std::string input_bam, 
		 // const std::string normal_bam,
		 const std::string prefix,
		 const int insert_size,
		 const bool scanall) 
{
  double walltime0 = get_wall_time();
  double cputime0 = get_cpu_time();

  std::cout << "[Step 1: Building BP-graph]\n";

  BamTools::RefVector references;
  
  std::unordered_map < std::string , int > id_ref_map;

  std::map < Position , std::map < Position , std::map < std::string , int > > >  BP_graph;

  //////////////////////////////////////////////////////////////////////
  //
  // Get break points
  //

  GetReferenceData ( input_bam, references, id_ref_map);

  //////////////////////////////////////////////////////////////////////
  //
  // Finding split reads related with BP candidates
  //
  
  Find_BP_paths(input_bam, id_ref_map, BP_graph, insert_size, scanall);

  /*
  //////////////////////////////////////////////////////////////////////
  //
  // Remove germline
  //
  if ( normal_bam.size() > 0 ) {
    remove_germline(normal_bam, id_ref_map, BP_graph, insert_size, scanall);
  }
  */

  //////////////////////////////////////////////////////////////////////
  //
  // Writing scores
  //

  Writing_files(prefix, BP_graph, references );

  //////////////////////////////////////////////////////////////////////
  //
  // Printing running time
  //

  double cputime1 = get_cpu_time();
  double walltime1 = get_wall_time();

  std::cout << "CPU time: " << cputime1 - cputime0 << " sec\n";
  std::cout << "Wall-clock time: " << walltime1 - walltime0 << "sec \n\n";

}


