//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#ifndef PEAK_MEMORY_MEARSURE
#define PEAK_MEMORY_MEARSURE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

/**
 * @brief Parsing line
 * @param [in] line
 */
std::size_t parseLine(char* line);

/**
 * @brief Peak memory use
 * @return Peak memory use in KB
 */
std::size_t get_peak_memory();

#endif
