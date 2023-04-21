#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <iostream>

int main ( int argc , char ** argv ){
  if ( argc != 4 ){
    std::cout << "Usage: etching_make_training_table etching_merge.vcf output_table.txt cutoff\n";
    return 1;
  }

  std::string infile ( argv[1] );
  std::string outfile ( argv[2] );
  int cutoff ( atoi ( argv[3] ) );

  std::ifstream fin ( infile.c_str() );

  if ( ! fin ){
    std::cerr << "ERROR!!! There is no input file: " << infile << "\n";
    std::cout << "Usage: etching_make_training_table etching_merge.vcf output_table.txt cutoff\n";
    return 1;
  }

  std::string tmp;
  std::string tmp1;

  std::ofstream fout ( outfile.c_str() );

  fout << "CR\tSR\tPE\tMQ\tDEPDIF\tTCB\tTorF\n";
  
  while ( std::getline ( fin , tmp ) ){
    if ( tmp[1] != '#' ) {
      break;
    }
  }

  std::string info_line;
  std::string feature;
  std::string detect;
  std::vector < std::string > ev(10);
  
  while ( std::getline ( fin , tmp ) ){
    std::stringstream ss ( tmp );
    
    for ( std::size_t i = 0 ; i < 7 ; i ++ ) std::getline ( ss , tmp1 , '\t' );

    std::getline ( ss , info_line , '\t' );

    std::getline ( ss , tmp1 , '\t' );
    std::getline ( ss , tmp1 , '\t' );
    
    if ( tmp1 != "NONE" ){
      std::stringstream ess ( tmp1 );
      std::string tmp2;
      
      for ( std::size_t i = 0 ; i < 10 ; i ++ ) {
	std::getline ( ess , tmp2 , ':' ) ;
	ev[i] = tmp2;
      }
      // svid = ev[0];
      feature = ev[1] + "\t" + ev[2] + "\t" + ev[3] + "\t" + ev[4] + "\t" + ev[5] + "\t" + ev[7];

      int count(0);
      while ( std::getline ( ss , tmp1 , '\t') ){
	if ( tmp1 != "NONE" ) count ++;
      }
      if ( count >= cutoff ) detect="1";
      else detect="0";

      fout << feature << "\t" << detect << "\n";
    }
  }

  fin.close();
  fout.close();

  return 0;
}
