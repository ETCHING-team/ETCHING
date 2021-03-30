#include "my_vcf.hpp"
#include <cmath>

double get_wall_time();
double get_cpu_time();
std::size_t get_peak_memory();

VCF_CLASS remove_germline ( VCF_CLASS container_t, VCF_CLASS container_n ){
  VCF_CLASS container;
  Position Pos1;
  Position Pos2;
  // std::string strand;

  Position Pos2_n;
  // std::string strand_n;

  for ( auto i : container_t ){
    bool check = 1 ;
    for ( int x = -10 ; x <= 10 ; x ++ ){
      Pos1 = i.first;
      Pos1.second += x;
      if ( container_n.find(Pos1) != container_n.end() ){
	check = 0 ;
	break;
      }
    }
    if ( check == 1 ){
      for ( auto vcf : i.second ){
	container.insert(vcf);
      }
    }
    else{
      for ( auto vcf : i.second ){
	Pos2.first = container_t.id_ref_map[vcf.chr2];
	// strand = vcf.strand;
	check = true ;

	for ( int x = -10 ; x <= 10 ; x ++ ){
	  Pos2.second = vcf.pos2 + x;

	  for ( auto vcf_n : container_n[Pos1] ){
	    Pos2_n.first = container_t.id_ref_map[vcf_n.chr2];
	    Pos2_n.second = vcf_n.pos2;
	    // strand_n = vcf_n.strand;
	    // if ( Pos2 == Pos2_n && strand == strand_n ){
	    if ( Pos2 == Pos2_n ){
	      check = false;
	      break;
	    }
	  }
	}
	if ( check == 1 ){
	  container.insert(vcf);
	}
      }
    }
  }

  return container;
}

/////////////////////////////////////////////////////////////////////////////////////
//
// Main
//


int main ( int argc , char ** argv ){
  if ( argc != 4 ){
    std::cout << "Usage:  somatic_filter  tumor.vcf  normal.vcf  input.bam >  output.vcf\n";
    return 0;
  }

  double cputime0 = get_cpu_time();
  double walltime0 = get_wall_time();

  std::string tmp;
  VCF vcf;
  VCF_CLASS container_t;
  VCF_CLASS container_n;
  VCF_CLASS container;

  std::string tumor_vcf = argv[1];
  std::string normal_vcf= argv[2];
  std::string input_bam = argv[3];


  /////////////////////////////////////////////////////////////////////////////
  //                                                                           
  // Bam initializing                                                          
  //

  BamTools::BamReader reader;
  BamTools::RefVector references;
  BamTools::BamAlignment al;

  reader.Open(input_bam);
  if ( !reader.OpenIndex(input_bam + ".bai") ) {
    reader.CreateIndex();
  }
  references = reader.GetReferenceData();


  /////////////////////////////////////////////////////////////////////////////
  //
  // Read vcf
  //

  for ( std::size_t i = 0 ; i < references.size() ; i ++ ){
    container_t.id_ref_map[references[i].RefName] = i ;
    container_t.ref_id_map[i] = references[i].RefName ;
  }
  container_t.id_ref_map["."]=-1;
  container_t.ref_id_map[-1]=".";

  container_n=container_t;

  /////////////////////////////////////////////////////////////////////////////////////
  //
  // Reading tumor vcf
  //

  std::cerr << "Reading " << tumor_vcf << "\n";
  // container_t.read_vcf_file(tumor_vcf);
  std::string header;
  std::ifstream fin ( tumor_vcf);
  while ( std::getline ( fin , tmp ) ){
    header += tmp + "\n";;
    if ( tmp[1] != '#' ){
      break;
    }
  }
  while ( std::getline ( fin , tmp ) ){
    vcf.parse_etching(tmp);
    container_t.insert(vcf);
  }
  fin.close();

  /////////////////////////////////////////////////////////////////////////////////////
  //
  // Reading normal vcf
  //

  std::cerr << "Reading " << normal_vcf << "\n";
  // container_n.read_vcf_file(normal_vcf);
  fin.open ( normal_vcf );

  while ( std::getline ( fin , tmp ) ){
    if ( tmp[1] != '#' ){
      break;
    }
  }
  while ( std::getline ( fin , tmp ) ){
    vcf.parse_etching(tmp);
    container_n.insert(vcf);
  }
  fin.close();

  
  std::map < int , std::map < int , int > > control_pos;

  for ( auto & i : container_n ){
    control_pos[i.first.first][i.first.second]=1;
    for ( auto & j : i.second ) {
      control_pos[container_n.id_ref_map[j.chr2]][j.pos2]=1;
    }
  }
  

  /////////////////////////////////////////////////////////////////////////////////////
  //
  // Remove germline
  //
  
  std::cerr << "Removing germline SVs";
  //container = remove_germline ( container_t, container_n );
  Position Pos1;
  Position Pos2;
  // std::string strand;

  Position Pos2_n;
  // std::string strand_n;

  for ( auto & i : container_t ){
    bool check = 1 ;
    for ( int x = -10 ; x <= 10 ; x ++ ){
      Pos1 = i.first;
      Pos1.second += x;

      if ( control_pos[Pos1.first].find(Pos1.second) != control_pos[Pos1.first].end() ){
	check = 0 ;
	// break;
      }
    }
    if ( check == 1 ){
      for ( auto & vcf : i.second ){
	container.insert(vcf);
      }
    }
    else{
      for ( auto & vcf : i.second ){
	bool check = true ;

	for ( int x = -10 ; x <= 10 ; x ++ ){
	  Pos2.first = container_n.id_ref_map[vcf.chr2];
	  Pos2.second = vcf.pos2 + x;

	  if ( control_pos[Pos2.first].find(Pos2.second) != control_pos[Pos2.first].end() ){
	    check = 0;
	    // break;
	  }
	}

	if ( check == 1 ){
	  container.insert(vcf);
	}
      }
    }
  }


  
  /////////////////////////////////////////////////////////////////////////////////////
  //
  // Print result
  //
  std::cout << header;
  container.write_no_header();


  //////////////////////////////////////////////////////////////////////////////////////////////
  //
  // End of program
  //

  double walltime1 = get_wall_time();
  double cputime1 = get_cpu_time();
  auto end = std::chrono::system_clock::now();
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);

  std::cerr << "===============================================================\n";
  std::cerr << "Total CPU time: " << cputime1 - cputime0 << " sec\n";
  std::cerr << "Total wall-clock time: " << walltime1 - walltime0 << " sec\n";
  std::cerr << "Peak memory: " << std::to_string((double)get_peak_memory()/(1024*1024)) << " GB\n";
  std::cerr << std::ctime(&end_time) << "\n";
  std::cerr << "===============================================================\n";



  return 0;
}
