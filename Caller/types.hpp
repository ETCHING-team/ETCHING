//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------


#ifndef MY_TYPES
#define MY_TYPES

#include <string>
#include <map>

using Position = std::pair < int , int > ;
using Feature = std::map < std::string , double > ;

typedef int64_t pos_type;
typedef std::tuple<pos_type,pos_type,pos_type> BP_type;
typedef std::pair<BP_type,BP_type> BND_type;

#endif

