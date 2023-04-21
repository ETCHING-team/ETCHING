#include "etching_info.hpp"
#include <sstream>
#include <string>

std::string get_program_name(int argc , char ** argv ){
  std::stringstream ss ( argv[0] );
  std::string token;
  while ( std::getline ( ss , token , '/' ));
  return token;
}
