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

#include "spp.h"

#include "Peak_memory_measure.hpp"
#include "CPU_TIME_measure.hpp"
#include <gzstream.h>



using namespace std;
using spp::sparse_hash_set;


#define TWO_BITS_MASK (3)
#define BITS_PER_BYTE (8)
#define BIG_ENOUGH (1024)


string decode(uint64_t encoded, bool rna_flag, size_t kmer_size);
uint64_t encode(const char *seq);
void encode_forw(uint64_t & value, const char & base, const size_t &kl);
uint64_t revcom(unsigned char x);
void encode_reve(uint64_t & value, const char & base, const size_t & kl);
uint64_t select_canonical (uint64_t & value1,uint64_t & value2);

uint64_t return_encode(string & kmer);
bool search_kmer(string & read1, string & qual1, string & read2, string & qual2);

void read_fastq();
void read_fastq_gz();
void encoding();
void write_fastq();
void SWAP();

int collector (string kmer_table, int kmer_size, string out_prefix, string read_file_1, string read_file_2, int gz_check, int NT);
   
#endif

