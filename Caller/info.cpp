#include "etching_caller.hpp"

void print_version(){
  std::cout << version << "\n";
}

void caller_usage(){
  std::cout 
    << "\n"
    << PROGRAM << "\n"
    << VERSION << "\n"
    << "\n"
    << "Usage:\tetching_caller [options] -b input.bam -o output_prefix -g reference.fa [options]\n"
    << "\n"
    << "Required:\n"
    << "\t" << "-b (string)\t" << "input bam file (required)\n"
    << "\t" << "-o (string)\t" << "Prefix of output file (required)\n"
    << "\t" << "-g (string)\t" << "Reference genome in fasta format (required)\n"
    << "\n"
    << "Procedure options:\n"
    // << "\t" << "-t (int)   \t" << "Number of threads [8]\n"
    << "\t" << "-B         \t" << "Print only BND without identification of SV types [NONE]\n"
    // << "\t" << "-S (int)   \t" << "Step: 1 from first step, 2 from second step, 3, from third step [1]\n"
    // << "\t" << "-i (int)   \t" << "Prefix of input file for -S 2 and -S 3\n"
    // << "\t" << "           \t" << "Do not use this option without -f option\n"
    << "\t" << "-A         \t" << "Using all split-reads in finding BPs [NONE]\n"
    // << "\t" << "-n (string)\t" << "Matched normal bam file for somatic call [NONE]\n"
    // << "\t" << "           \t" << "Note!!! Do not need to use this option, if you used ETCHING-FILTER\n"
    << "\n"
    << "Information of sequencing data:\n"
    // << "\t" << "-D (float) \t" << "Sequencing depth [50 for WGS, 100 for WES, 500 for PANEL]\n"
    << "\t" << "-D (float) \t" << "Sequencing depth [50]\n"
    << "\t" << "-P (float) \t" << "Tumor purity [0.75]\n"
    // << "\t" << "-T (string)\t" << "Type of data: One of WGS, WTS, WES, PANEL, and so on [WGS]\n"
    << "\t" << "-I (int)   \t" << "Insert-size [500]\n"
    << "\t" << "-O (string)\t" << "Read orientation: FR for Paired-End, and RF for Mate-Pair [FR]\n"
    << "\n"
    << "Information of this program:\n"
    << "\t" << "-h         \t" << "Print help\n"
    << "\n"
    << CONTACT << "\n"
    << "\n";
};
