//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#ifndef MAKE_GENE_LIST
#define MAKE_GENE_LIST

#include <iostream>
#include <fstream>
#include <string>
#include <map>

void parsing_gff_info(std::string input, std::map < std::string , std::string > & info );
int64_t make_gene_list ( std::string arg1, std::string arg2, std::string arg3, std::string arg4);

#endif
