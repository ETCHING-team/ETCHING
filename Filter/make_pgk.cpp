#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>

void usage(){
  std::cout << "Usage: make_pgk [options] -i reference_fasta.list [options]\n"
	    << "\n"
	    << "\t" << "-o <string>\t" << "Output prefix [PGK]\n"
	    << "\t" << "-v <string>\t" << "SNV/indel vcf for control. For instace, dbSNP.vcf\n"
	    << "\t" << "-g <string>\t" << "Reference genome for input vcf\n"
	    << "\t" << "-l <int>   \t" << "k-mer length [31]\n"
	    << "\t" << "-D <int>   \t" << "Directory of KMC-3 [null]\n"
	    << "\t" << "-t <int>   \t" << "Number of threads [8]\n";
}

void printf_alt_kmer ( std::map < std::string , std::string > & genome , std::string & id , int & pos , std::string & alt , std::string & variant_id, std::string & REF, int & kl , int & count, std::ofstream & fout){
  std::string ID=std::to_string(kl) + "-mer_";
  std::stringstream ss ( alt );
  std::string ref;
  std::string tmp;
  std::vector < std::string > alt_vec;
  while ( std::getline ( ss , tmp , ',' ) ){
    alt_vec.push_back (tmp);
  }

  int start = pos - kl;
  int end = pos + kl + REF.size() - 1;
  int shift = 0 ;
  if ( start < 0 ){
    shift = start;
    start = 0 ;
  }
  if ( end > (int) genome[id].size() ){
    end = (int) genome[id].size() ;
  }

  ref = genome[id].substr ( start, end - start - 1 );
  
  for ( auto ALT : alt_vec ){
    std::string output;
    output = ref.substr ( 0 , kl - 1 + shift ) + ALT + ref.substr(kl + REF.size() -1 );
    fout << ">" << ID << count << " " << id << " " << pos << "\t" << REF << "\t" << ALT <<  "\n" << output << "\n";
    count ++;
  }
}





void extract_seq_related_with_snv_indel( std::string genome_file , std::string control_vcf , std::string snvindel, int kl){

  std::map < std::string , std::string > genome;

  // reading genome file
  std::cerr << "\tReading genome: " + genome_file + "\n";

  std::string id;
  std::string seq;
  std::string tmp;
  
  std::ifstream fin ( genome_file.c_str() );
    
  std::getline ( fin , tmp );
  tmp = tmp.substr(1);
  std::stringstream ss0 ( tmp );
  ss0 >> id;

  while ( std::getline ( fin , tmp ) ){
    if ( tmp[0] != '>' ){
      seq += tmp;
    }
    else {
      genome[id]=seq;
      seq.clear();
	
      tmp = tmp.substr(1);
      std::stringstream ss ( tmp );
      ss >> id;

    }
  }
  genome[id]=seq;
  seq.clear();
    
  fin.close();
    
  // reading vcf file and extract sequence associated with snv or indel
    
  std::cerr << "\tExtract sequences : " + control_vcf + "\n";

  int pos;
  std::string alt;
  std::string variant_id;
  std::string REF;
  int count = 1;

  std::ofstream fout ( snvindel.c_str() );
    
  fin.open ( control_vcf.c_str() );
  while ( std::getline ( fin , tmp )){
    if ( tmp.substr(0,2) != "##" ){
      break;
    }
  }
  while ( std::getline ( fin , tmp ) ){
    std::stringstream ss ( tmp );
    ss >> id >> pos >> variant_id >> REF >> alt;
    if ( genome.find(id) != genome.end() ){
      printf_alt_kmer ( genome , id , pos , alt , variant_id, REF, kl , count, fout);
    }
  }
  fin.close();
    
  fout.close();
}


int main ( int argc, char ** argv ){

  if ( argc == 1 ){
    usage();
    return 0;
  }

  std::string fasta_list;
  std::string prefix="PGK";

  std::string control_vcf;
  std::string genome_file;

  int kl = 31;

  std::string DIR;
  int num_threads=8;
  int maxRAM=12;

  int opt;

  while ( (opt = getopt ( argc, argv, "i:o:v:g:l:D:t:m:h" ) ) != -1 ){

    switch ( opt ) {
      
    case 'i': fasta_list = optarg; break;
    case 'o': prefix = optarg; break;
      
    case 'v': control_vcf = optarg; break;
    case 'g': genome_file = optarg; break;
      
    case 'l': kl = atoi(optarg); break;

    case 'D': DIR = optarg ; break;
    case 't': num_threads = atoi ( optarg ); break;
    case 'm': maxRAM = atoi ( optarg ); break;

    case 'h': usage(); return 0 ;
      
    default: std::cout << "\tInvalid option\n"; usage(); return 1;
      
    }
  }
  
  // check options
  if ( fasta_list.size() == 0 ){
    std::cout << "ERROR!!! -i option is required.\n";
    return -1;
  }

  
  // pop_back if DIR has "/" charactor at the tail
  if ( DIR.size() != 0 ){
    if ( DIR[DIR.size()-1] != '/' ){
      DIR += '/';
    }
  }

  // check kmc
  std::string kmc_path=DIR+"kmc > /dev/null";
  if ( system(kmc_path.c_str()) > 0 ){
    std::cout << "ERROR!!! We cannot find kmc. Please set -D option for kmc, or install KMC-3.";
    return -1;
  }

  return 0;

  std::string tmp;

  // update genome_vec for storing fasta files

  std::vector < std::string > genome_vec;
  std::ifstream fin ( fasta_list.c_str() );
  while ( fin >> tmp ){
    genome_vec.push_back(tmp);
  }
  fin.close();

  std::cerr << "[Controls for filtration]\n";
  for ( auto i : genome_vec ){
    std::cerr << "\t" << i << "\n";
  }
  std::cerr << "\n";

  // if contorl_vcf was given, extract sequences related with the vcf
  if ( control_vcf.size() * genome_file.size() > 0 ){
    std::cerr << "[Extract sequences related with SNVs or INDELs]\n";

    std::string snvindel = control_vcf + ".SNV_INDEL.fasta";
    genome_vec.push_back(snvindel);
    extract_seq_related_with_snv_indel ( genome_file , control_vcf , snvindel , kl );
    std::cerr << "\n";
  }


  // build kmc database
  
  std::string updated_fasta_list_file = prefix + ".reference.list";
  std::ofstream fout ( updated_fasta_list_file.c_str() );
  for ( auto i : genome_vec ){
    fout << i << "\n";
  }
  fout.close();

  std::string command ;
  std::string tmp_dir = prefix + ".tmp_dir";
  command = "mkdir -p " + tmp_dir;
  system(command.c_str());

  std::cerr << "[Build k-mer database]\n";
  
  std::string output_log = prefix + ".log";
  command
    = DIR + "kmc -k" + std::to_string(kl) + " -t" + std::to_string(num_threads) 
    + " -m" + std::to_string ( maxRAM ) + " -v -ci1 -fm @" + updated_fasta_list_file + " " + prefix + " " + tmp_dir + " > " + output_log + " 2>&1" ;
  system(command.c_str());

  command = "rm -rf " + tmp_dir;
  system(command.c_str());


  std::cerr << "[Output]\n";
  std::cerr << "\t" << prefix << ".kmc_pre\n";
  std::cerr << "\t" << prefix << ".kmc_suf\n";

  return 0;
}











