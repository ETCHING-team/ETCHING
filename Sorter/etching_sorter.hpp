//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#ifndef ETCHING_SORTER
#define ETCHING_SORTER

#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <unistd.h>
#include "types.hpp"

Feature harmonic_average (Feature f1, Feature f2);
void print_header(std::ofstream & fout);
void print_combi_feat(std::string id, std::string type, Feature f, std::vector < std::string > feat_vec, std::ofstream & fout);
std::string return_svtype(std::string info);
std::string cut_str(std::string & input, std::string key);
Feature parse_feature ( std::string format , std::string format1);
void calc_feature ( std::string infile, std::string outfile );
void razor ( std::string infile, std::string score_file, std::string outfile, double cutoff, std::string method, int tagging);

#endif
