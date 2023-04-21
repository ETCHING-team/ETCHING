//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#ifndef CPU_TIME_MEASURE
#define CPU_TIME_MEASURE

#include <time.h>
#include <sys/time.h>

/**
 * @brief Get wall-clock time to calculate running time (wall-clock)
 * @return double precision flaot
 */
double get_wall_time();

/**
 * @brief Get CPU time to calculate running time (CPU)
 * @return double precision flaot
 */
double get_cpu_time();

#endif

