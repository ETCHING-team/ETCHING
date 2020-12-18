//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#include "etching_caller.hpp"
#include "my_vcf.hpp"


void bp_to_vcf ( const std::string genome,
                 const std::string input_bam,
                 const std::string infile_pref,
                 const std::string prefix,
                 const int insert_size,
                 const bool typing,
                 const std::string read_orientation,
                 const double tumor_purity,
                 const double sequencing_coverage,
                 const std::string data_type)
{
  std::cout << "[Step 2: Printing BPs]\n";

  BamTools::BamReader reader;
  BamTools::RefVector references;

  reader.Open(input_bam);
  if ( !reader.OpenIndex(input_bam + ".bai") ) {
    reader.CreateIndex();
  }
  references = reader.GetReferenceData();

  
  
  BamTools::BamAlignment al;
  int read_length=0, count = 0, max = 10000;
  while ( count ++ < max ){
    reader.GetNextAlignment(al);
    read_length < al.Length ? read_length = al.Length : read_length = read_length ;
  }
  
  VCF_CLASS container;

  
  container.refvector=references;

  
  for ( std::size_t i = 0 ; i < references.size() ; i ++ ){
    container.id_ref_map[references[i].RefName] = i ;
    container.ref_id_map[i] = references[i].RefName ;
  }
  container.id_ref_map["."]=-1;
  container.ref_id_map[-1]=".";

  const std::string pair_input = infile_pref + ".pair.txt" ;

  double cputime0 = get_cpu_time(), walltime0 = get_wall_time();
  std::cout << "Reading reference genome: " << genome << "\n";
  container.get_genome(genome);


  /////////////////////////////////////////////////////////////////////////////////////
  //
  // Reading Clipped read count
  //

  std::cout << "Searching candidates\n";

  std::string chr1, chr2, strand, direction;
  int pos1, pos2, tmp_int;
  
  std::ifstream fin;
  
  fin.open ( pair_input.c_str() );
  count = 1;
  while ( fin >> chr1 >> pos1 >> chr2 >> pos2 >> strand >> tmp_int ){
    
    Position Pos1, Pos2;
  
    Pos1.first = container.id_ref_map[chr1];
    Pos1.second= pos1;
    
    Pos2.first = container.id_ref_map[chr2];
    Pos2.second= pos2;
    
    std::string strand_mate;

    if ( strand == "FR" ) strand_mate = "RF";
    else if ( strand == "RF" ) strand_mate = "FR";
    else if ( strand == "FF" ) strand_mate = "FF";
    else if ( strand == "RR" ) strand_mate = "RR";
    else strand_mate = ".";

    if ( container.check_vcf( chr1, pos1, chr2, pos2, strand) == 0 && ( chr1 != chr2 || pos1 != pos2 ) ){
      if ( chr2 != "." ){
	std::string sv_id   = "BND" + std::to_string(count++);
	std::string mate_id = "BND" + std::to_string(count++);
	container.insert ( chr1, pos1, chr2, pos2, sv_id, mate_id, strand, tmp_int, "BND");
	container.insert ( chr2, pos2, chr1, pos1, mate_id, sv_id, strand_mate, tmp_int, "BND");
      }
      else{
	std::string sv_id   = "SND" + std::to_string(count++);
	container.insert ( chr1, pos1, chr2, pos2, sv_id, ".", strand, tmp_int, "SND");
      }
    }
  }
  
  fin.close();

  container.genome.clear();

  // global features
  for ( auto & i : container ){
    for ( auto & j : i.second ){
      j.seqdep = sequencing_coverage;
      j.purity = tumor_purity;
      j.data_type = data_type;
    }
  }
  
  /////////////////////////////////////////////////////////////////////////////////////
  //
  // calc_feature; 
  //
  
  std::cout << "Calculating features\n";

  // input local features
  container.calc_features ( input_bam, read_length, insert_size, 1);
  // container.calc_features_general ( input_bam, read_length, insert_size, 1);
  // container.fill_mate_feature ();

  /////////////////////////////////////////////////////////////////////////////////////
  //
  // Typing SV
  //

  container.make_header();
  container.make_info();

  VCF_CLASS container_SV;
  if ( typing) container_SV=typing_SV(container);
  
  /////////////////////////////////////////////////////////////////////////////////////
  //
  // Writing vcf files
  //

  std::string outfile_vcf = prefix + ".BND.vcf";
  
  container.fwrite(outfile_vcf);

  if ( typing ){
    outfile_vcf = prefix + ".SV.vcf";
    container_SV.fwrite(outfile_vcf);
  }
  
  double cputime1 = get_cpu_time();
  double walltime1 = get_wall_time();
  std::cout << "CPU time: " << cputime1 - cputime0 << " sec\n";
  std::cout << "Wall-clock time: " << walltime1 - walltime0 << " sec\n\n";

}

