//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#ifndef FIND_FUSION_GENE
#define FIND_FUSION_GENE

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>


using Position = std::pair < std::string , int > ;
using Loci = std::pair < int , int > ;

void read_bp_pair(const std::string BP_pair_file, std::map < int , std::map < Position , std::map < Position , std::string > > > & BP_pair);

void read_gene_list(const std::string gene_list_file,
		    std::map < std::string , std::map < Loci , std::string > > & gene_list_for,
		    std::map < std::string , std::map < Loci , std::string > > & gene_list_rev);

void find_fusion_gene ( std::map < int , std::map < Position , std::map < Position , std::string > > > & BP_pair,
			std::map < std::string , std::map < Loci , std::string > > & gene_list_for,
			std::map < std::string , std::map < Loci , std::string > > & gene_list_rev,
			std::set < std::string > & fusion_genes, int window, int strand_aware );

int find_fg ( std::string BP_pair_file, std::string gene_list_file, int window, int strand_aware );

#endif
