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
using FEATURE = std::map < std::string , Feature > ;

#endif

