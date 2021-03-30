#ifndef LOCAL_MIN_KMER_DEPTH
#define LOCAL_MIN_KMER_DEPTH

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

int main ( int argc , char ** argv ){
  
  std::vector <std::size_t> depth;
  std::vector <std::size_t> freq;
  std::size_t d,f;
  std::size_t local_min;
  
  std::ifstream fin(argv[1]);
  while ( fin >> d >> f ){
    depth.push_back(d);
    freq.push_back(f);
  }
  fin.close();
  
  for ( size_t i = 1 ; i < depth.size() ; i ++ ){
    if ( freq[i-1] < freq[i] ){
      local_min = depth[i-1];
      break;
    }
    // }
    // else if ( data_type == "P" ) { // for panel
    //   for ( size_t i = 1 ; i < depth.size() ; i ++ ){
    // 	if ( freq[0] > 10 * freq[i] ){
    // 	  // local_min = depth[i];
    // 	  if ( freq[i-1] - freq[0]/10 < freq[0]/10 - freq[i] ){
    // 	    local_min = depth[i-1];
    // 	  }
    // 	  else{
    // 	    local_min = depth[i];
    // 	  }
    // 	  break;
    // 	}
    //   }
    if ( local_min < 2 ) local_min = 2 ;
  }
  std::cout << local_min << "\n";

  return 0;
}
#endif
