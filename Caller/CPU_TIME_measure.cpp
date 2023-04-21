//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#include "CPU_TIME_measure.hpp"

double get_wall_time(){
  struct timeval time;
  if (gettimeofday(&time,NULL)){

    return 0;
  }
  return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

double get_cpu_time(){
  return (double)clock() / CLOCKS_PER_SEC;
}

