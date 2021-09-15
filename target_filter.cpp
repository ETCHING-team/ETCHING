#include "my_vcf.hpp"
#include <sstream>

double get_wall_time();
double get_cpu_time();
std::size_t get_peak_memory();

inline bool is_in_target ( std::string line , std::map < std::string , std::map < int , int > > & target_region ){
  std::stringstream ss ( line );
  std::string chr1, chr2, tmp, alt, info;
  int pos1(-1), pos2(-1);
  ss >> chr1 >> pos1 >> tmp >> tmp >> alt >> tmp >> tmp >> info;
  if ( target_region.find ( chr1 ) != target_region.end() ){
    for ( auto i : target_region[chr1] ){
      if ( pos1 >= i.first && pos1 <= i.second ){
	return 1;
      }
    }
  }

  // search for chr2:pos2

  std::stringstream info_ss ( info );
  std::map < std::string , std::string > info_map;
  std::string key;
  std::string val;
  
  while ( std::getline ( info_ss , tmp , ';' ) ){
    std::stringstream tmp_ss ( tmp );
    std::getline ( tmp_ss , key , '=' );
    std::getline ( tmp_ss , val , '=' );
    info_map[key]=val;
  }
  
  if ( info_map.find("END") != info_map.end () ){
    pos2=atoi(info_map["END"].c_str());
    if ( info_map.find("CHR2") != info_map.end () ){
      chr2=info_map["CHR2"];
    }
    else{
      chr2=chr1;
    }
  }
  else{
    std::size_t found1 = ( alt.find ( ":" ) != std::string::npos ) ;
    std::size_t found2 = ( alt.find ( "[" ) != std::string::npos ) ;
    std::size_t found3 = ( alt.find ( "]" ) != std::string::npos ) ;

    if ( found1 && found2 && ! found3 ){
      std::stringstream ss1(alt);
      std::getline ( ss1 , tmp , '[' );      
      std::getline ( ss1 , alt , '[' );      
    }
    else if ( found1 && ! found2 && found3 ){
      std::stringstream ss1(alt);
      std::getline ( ss1 , tmp , ']' );      
      std::getline ( ss1 , alt , ']' );      
    }
    else return 0;

    std::stringstream ss1(alt);
    std::getline ( ss1 , chr2 , ':' );
    std::getline ( ss1 , tmp , ':' );
    pos2=atoi(tmp.c_str());
  }


  //////////////
  /*
  std::size_t found = line.find ( "CHR2=" );
  if ( found != std::string::npos ){
    line = line.substr ( found + 5 );
    found = line.find ( ";" );
    chr2 = line.substr ( 0, found );
  
    found = line.find ( "END=" );
    if ( found != std::string::npos ){
      line = line.substr ( found + 4 );
      found = line.find ( ";" );
      pos2 = atoi( line.substr ( 0, found ).c_str() );  
    }
  }
  else {
    found = alt.find ( ":" );
    if ( found != std::string::npos ) {
      if ( alt[0]==']' || alt[0]=='[' ){
	alt = alt.substr(1,alt.size()-3);
      }
      else if ( alt[1]==']' || alt[1]=='[' ) {
	alt = alt.substr(2,alt.size()-3);	
      }
      else return 0;
    }
    else return 0;
    std::stringstream ss1(alt);
    std::getline ( ss1 , chr2 , ':' );
    std::getline ( ss1 , tmp , ':' );
    pos2=atoi(tmp.c_str());
  }
  */

  if ( target_region.find ( chr2 ) != target_region.end() ){
    for ( auto i : target_region[chr2] ){
      if ( pos2 >= i.first && pos2 <= i.second ){
	return 1;
      }
    }
  }

  return 0;
}

void target_filter_usage(){
  std::cerr << "Usage: target_filter input.vcf target.bed > output.vcf\n";
}

int main ( int argc , char ** argv ){
  if ( argc != 3 ){
    target_filter_usage();
    return 1;
  }
  
  std::string infile_vcf(argv[1]);
  std::string infile_bed(argv[2]);

  std::string tmp;
  tmp = "[Reading bed: " + infile_bed + "]\n";
  std::cerr << tmp ;
  std::ifstream fin ( infile_bed.c_str() );
  std::string chr;
  int pos1;
  int pos2;
  std::map < std::string , std::map < int , int > > target_region;
  while ( std::getline ( fin , tmp ) ){
    if ( tmp[0] != '#' ){
      std::stringstream ss ( tmp );
      ss >> chr >> pos1 >> pos2;
      target_region[chr][pos1]=pos2;
    }
  }
  fin.close();

  std::cerr << "[Reading vcf: " << infile_vcf << "]\n";
  fin.open ( infile_vcf );
  while ( std::getline ( fin , tmp ) ){
    std::cout << tmp << "\n";
    if ( tmp[1] != '#' ){
      break;
    }
  }
  while ( std::getline ( fin , tmp ) ){
    if ( is_in_target ( tmp , target_region ) ) {
      std::cout << tmp << "\n";
    }
  }
  fin.close();

  std::cerr << "[Finished]\n" ;

  return 0;
}
