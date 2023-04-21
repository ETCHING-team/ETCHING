#include "etching_caller.hpp"

void caller_usage(int argc , char ** argv){
  std::cout 
    << "\n"
    << get_program_name(argc, argv) << "\n"
    << "Version: " << ETCHING_VERSION << " (" << RELEASE_DATE << ")\n"
    << "\n"
    << "Usage:\tetching_caller [options] -b input.bam -o output_prefix -g reference.fa [options]\n"
    << "\n"
    << "Required:\n"
    << "\t" << "-b (string)\t" << "input bam file (required)\n"
    << "\t" << "-o (string)\t" << "Prefix of output file (required)\n"
    << "\t" << "-g (string)\t" << "Reference genome in fasta format (required)\n"
    << "\n"
    << "Options:\n"
    << "\t" << "-I (int)   \t" << "Insert-size [500]\n"
    << "\t" << "-O (string)\t" << "Read orientation: FR for Paired-End, and RF for Mate-Pair [FR]\n"
    << "\t" << "-B         \t" << "Write only BND without identification of SV types\n"
    << "\t" << "-R         \t" << "Rescue SVs not supported by split-reads\n"
    << "\t" << "--bp-pair  \t" << "Write only BP-pairs\n"
    << "\t" << "-h         \t" << "Print help\n"
    << "\n"
    << ETCHING_CONTACT << "\n"
    << "\n";
};

