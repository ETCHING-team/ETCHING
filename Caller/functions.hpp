//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------


#ifndef MY_FUNCTIONS
#define MY_FUNCTIONS

#include <string>
#include <vector>
#include <unordered_map>
#include <api/BamAlignment.h>
#include "types.hpp"


int check_concord ( const Position & Pos1, const Position & Pos2, const int & insert_size );
int FindConcordant ( const Position & Pos_read , const Position & Pos_split , const Position & Pos_mate , const int & insert_size );
std::string parse_SA_tag ( std::string & Tag );
Position find_split ( const std::string & tag , const std::unordered_map < std::string , int > & id_ref_map );
std::string mate_strand ( const std::string & tag );
std::vector < BamTools::CigarOp > cigar_parse ( std::string cigar );
std::string currentDate();
std::string currentDateTime();


#endif

