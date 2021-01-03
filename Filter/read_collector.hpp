//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#ifndef READ_COLLECTOR
#define READ_COLLECTOR

#include <unistd.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <omp.h>
#include <vector>
#include <utility>
#include <thread>
#include <bitset>

#include "Peak_memory_measure.hpp"
#include "CPU_TIME_measure.hpp"

#include <gzstream.h>

using namespace std;

int64_t collector (string kmer_table, int64_t kmer_size, string out_prefix, string read_file_1, string read_file_2, int64_t gz_check, int64_t NT);
int64_t collector_id_mode_single (string kmer_table, int64_t kmer_size, string out_prefix, string read_file_1, int64_t gz_check, int64_t NT);
int64_t collector_id_mode (string kmer_table, int64_t kmer_size, string out_prefix, string read_file_1, string read_file_2, int64_t gz_check, int64_t NT);
int64_t collector_bam (string kmer_table, int64_t kmer_size, string out_prefix, string read_file_b, int64_t NT);
int64_t collector_single (string kmer_table, int64_t kmer_size, string out_prefix, string read_file_1, int64_t gz_check, int64_t NT);
   
#endif

