//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#include "Peak_memory_measure.hpp"

/**
 * @brief Parsing a line to get memory use
 * @return Memory use in KB
 */

inline std::size_t parseLine(char* line){
  std::size_t i = strlen(line);
  while (*line < '0' || *line > '9') line++;
  line[i-3] = '\0';
  i = atoi(line);
  return i;
}

std::size_t get_peak_memory()
{
  FILE* file = fopen("/proc/self/status", "r");
  std::size_t result = 0;
  char line[128];
  
  while (fgets(line, 128, file) != NULL)
    {
      if (strncmp(line, "VmPeak:", 6) == 0)
        {
	  result = parseLine(line);
	  break;
        }
    }
  fclose(file);
  return result;
}


