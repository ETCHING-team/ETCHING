//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#ifndef ETCHNING_CUTTER_FUNCTION
#define ETCHNING_CUTTER_FUNCTION

#include "etching_sorter.hpp"


Feature harmonic_average (Feature f1, Feature f2){
  for ( auto & i : f1 ){
    if ( i.second * f2[i.first] != 0 ){
      i.second = ( 2 * i.second * f2[i.first] ) / ( i.second + f2[i.first]) ;
    }
    else i.second = 0;
  }
  return f1;
}


std::string cut_str(std::string & input, std::string key){
  std::size_t found = input.find(key);
  std::string output;
  if ( found != std::string::npos ){
    output = input.substr(0,found);
    input = input.substr(found+key.size());
  }
  else{
    output = input;
    input.clear();
  }
  return output;
}                                    

std::string return_svtype(std::string info){
  std::string key = "SVTYPE=";
  std::size_t found = info.find(key);
  return info.substr(found+key.size(),3);
} 

Feature parse_feature ( std::string format , std::string format1){
  Feature output;

  std::string name;
  std::string tmp;
  double value;
  std::string key=":";
  std::size_t sz;

  while ( format.size() != 0 ){
    name = cut_str(format,key);
    tmp = cut_str(format1,key);
    value = std::stod(tmp,&sz);
    output[name]=value;
  }

  return output;
}


#endif
