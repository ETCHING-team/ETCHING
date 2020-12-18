//--------------------------------------------------------------------
// Copyright 2018
// Written by Jang-il Sohn (sohnjangil@gmail.com)
// Bioinformatic and Genomics Lab., Hanyang University, Seoul, Korea
// Released under the MIT license
//--------------------------------------------------------------------

#include "read_collector.hpp"

#include "encode.hpp"

class container{
private:
public:
  container(){};
  container(int l){
    init(l);
  };
  ~container(){
    clear();
  };
  int num_threads;
  vector<string> id_vec1;
  vector<string> read_vec1;
  vector<string> qual_vec1; 
  vector<string> id_vec2;
  vector<string> read_vec2;
  vector<string> qual_vec2;
  vector<bool> check_filter;
  size_t count_read;

  void init(int l){
    id_vec1.resize(l);
    read_vec1.resize(l);
    qual_vec1.resize(l);

    id_vec2.resize(l);
    read_vec2.resize(l);
    qual_vec2.resize(l);

    check_filter.resize(l); 

    count_read=0;
  }

  void clear(){
    id_vec1.clear();
    read_vec1.clear();
    qual_vec1.clear();
    id_vec2.clear();
    read_vec2.clear();
    qual_vec2.clear();
    check_filter.clear();
  }
};
container revolver[4];

  
spp::sparse_hash_set<uint64_t> kmerDB_set;

//size_t kmer_size;
size_t non_ref_kmer_cut(1);
//int num_threads;

int reading_chunk;
int chunk_size;

FILE* readFile_1;
FILE* readFile_2;

igzstream readFile_1_gz;
igzstream readFile_2_gz;

FILE* outfile_1_out;
FILE* outfile_2_out;

double read_wall_clock;
double read_cpu_time;

double encoding_wall_clock;
double encoding_cpu_time;

double write_wall_clock;
double write_cpu_time;



uint64_t return_encode(string & kmer){
  uint64_t forward(0);
  uint64_t reverse(0);
  size_t size=kmer.size();
  for ( size_t i = 0 ; i < size ; i ++ ){
    encode_forw(forward,kmer[i],size);
    encode_reve(reverse,kmer[i],size);
  }

  return forward > reverse ? forward : reverse;
}




bool search_kmer(string & read1, string & qual1, string & read2, string & qual2, size_t kmer_size)
{
  size_t checker1(0);
  size_t checker2(0);

  uint64_t encode_k(0);
  uint64_t reversed(0);

  uint64_t canonical;
  
  size_t Ncount(0) ;
  size_t i;

  /////////////////////////////////////////////////////////////////////////
  //
  // count non-ref. k-mers in read_1
  //
  
  if ( read1.size() >= kmer_size ){

    for ( i = 0 ; i < kmer_size ; i ++ ){
      encode_forw(encode_k,read1[i],kmer_size);
      encode_reve(reversed,read1[i],kmer_size);
      Ncount += read1[i] == 'N' ;
    }
    
    if ( Ncount == 0 ) {
      canonical = select_canonical( encode_k, reversed ) ;
      if ( kmerDB_set.find(canonical) != kmerDB_set.end() ) checker1 ++ ;
    }
    
    for ( ; i < read1.length() ; i++) {
      encode_forw(encode_k,read1[i],kmer_size);
      encode_reve(reversed,read1[i],kmer_size);
      Ncount += read1 [ i ] == 'N' ;
      Ncount -= read1 [ i - kmer_size ] == 'N' ;
      
      if ( Ncount == 0 ) {
	canonical = select_canonical( encode_k, reversed ) ;
	if ( kmerDB_set.find(canonical) != kmerDB_set.end() ) checker1 ++ ;
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////
  //
  // count non-ref. k-mers in read_2
  //

  if ( read2.size() >= kmer_size ){
    
    // must be initialized prior to check read_2
    encode_k = 0;
    reversed = 0;
    Ncount = 0  ;
    
    for ( i = 0 ; i < kmer_size ; i ++ ){
      encode_forw(encode_k,read2[i],kmer_size);
      encode_reve(reversed,read2[i],kmer_size);
      Ncount += read2[i] == 'N' ;
    }
    
    if ( Ncount == 0 ) {
      canonical = select_canonical( encode_k, reversed ) ;
      if ( kmerDB_set.find(canonical) != kmerDB_set.end() ) checker2 ++ ;
    }
    
    for ( ; i < read2.length() ; i++) {
      encode_forw(encode_k,read2[i],kmer_size);
      encode_reve(reversed,read2[i],kmer_size);
      Ncount += read2 [ i ] == 'N' ;
      Ncount -= read2 [ i - kmer_size ] == 'N' ;
      
      if ( Ncount == 0 ) {
	canonical = select_canonical( encode_k, reversed ) ;
	if ( kmerDB_set.find(canonical) != kmerDB_set.end() ) checker2 ++ ;
      }
    }
  }

  // return result
  return ( checker1 >= non_ref_kmer_cut || checker2 >= non_ref_kmer_cut ) ;
}



void find_filtered_read(size_t start, size_t end, int kmer_size)
{
  string read1;
  string qual1;
  string read2;
  string qual2;

  for ( size_t i = start ; i < end ; i ++ ){
    read1=revolver[1].read_vec1[i];
    qual1=revolver[1].qual_vec1[i];
    read2=revolver[1].read_vec2[i];
    qual2=revolver[1].qual_vec2[i];
    revolver[1].check_filter[i]=search_kmer(read1, qual1, read2, qual2, kmer_size);
  }
}



// void read_fastq(vector<string>  id_vec1, vector<string>  read_vec1, vector<string>  qual_vec1, 
// 		vector<string>  id_vec2, vector<string>  read_vec2, vector<string>  qual_vec2,
// 		vector<bool>  check_filter, size_t  count_read){
  
void read_fastq(){
  
  const double walltime0(get_wall_time());
  const double cputime0(get_cpu_time());

  size_t count_read = 0;

  char* id1 = NULL;
  char* read1 = NULL;
  char* desc1 = NULL;
  char* qual1 = NULL;
  char* id2 = NULL;
  char* read2 = NULL;
  char* desc2 = NULL;
  char* qual2 = NULL;

  size_t len11, len21;
  size_t len12, len22;
  size_t len13, len23;
  size_t len14, len24;

  while ( getline ( &id1 , &len11, readFile_1 ) != -1 && getline ( &id2 , &len21, readFile_2 ) != -1 &&
	  getline ( &read1 , &len12, readFile_1 ) != -1 && getline ( &read2 , &len22, readFile_2 ) != -1 &&
	  getline ( &desc1 , &len13, readFile_1 ) != -1 && getline ( &desc2 , &len23, readFile_2 ) != -1 &&
	  getline ( &qual1 , &len14, readFile_1 ) != -1 && getline ( &qual2 , &len24, readFile_2 ) != -1 )
    {

      revolver[0].id_vec1[count_read] = id1;
      revolver[0].read_vec1[count_read] = read1;
      revolver[0].qual_vec1[count_read] = qual1;

      revolver[0].id_vec2[count_read] = id2;
      revolver[0].read_vec2[count_read] = read2;
      revolver[0].qual_vec2[count_read] = qual2;

      // remove '\n'
      revolver[0].id_vec1[count_read].pop_back();
      revolver[0].read_vec1[count_read].pop_back();
      revolver[0].qual_vec1[count_read].pop_back();
      
      revolver[0].id_vec2[count_read].pop_back();
      revolver[0].read_vec2[count_read].pop_back();
      revolver[0].qual_vec2[count_read].pop_back();

      revolver[0].check_filter[count_read] = 0 ;

      count_read ++;
      
      if ( ! (count_read % reading_chunk ) ) {
	break;
      }
    }

  // std::cout << revolver[0].read_vec1[0] << "\n";
  revolver[0].count_read = count_read;

  if(id1) free(id1);
  if(id2) free(id2);
  if(read1) free(read1);
  if(read2) free(read2);
  if(desc1) free(desc1);
  if(desc2) free(desc2);
  if(qual1) free(qual1);
  if(qual2) free(qual2);

  const double walltime1(get_wall_time());
  const double cputime1(get_cpu_time());
  
  read_wall_clock += walltime1 - walltime0;
  read_cpu_time += cputime1 - cputime0;
}


void read_fastq_gz(){
  
  const double walltime0(get_wall_time());
  const double cputime0(get_cpu_time());
  
  size_t count_read = 0;

  string id1;
  string read1;
  string desc1;
  string qual1;
  string id2;
  string read2;
  string desc2;
  string qual2;

  while ( getline ( readFile_1_gz, id1 ) && getline ( readFile_2_gz, id2 ) &&
	  getline ( readFile_1_gz, read1 ) && getline ( readFile_2_gz, read2 ) &&
	  getline ( readFile_1_gz, desc1 ) && getline ( readFile_2_gz, desc2 ) &&
	  getline ( readFile_1_gz, qual1 ) && getline ( readFile_2_gz, qual2 ) )
    {

      revolver[0].id_vec1[count_read] = id1;
      revolver[0].read_vec1[count_read] = read1;
      revolver[0].qual_vec1[count_read] = qual1;

      revolver[0].id_vec2[count_read] = id2;
      revolver[0].read_vec2[count_read] = read2;
      revolver[0].qual_vec2[count_read] = qual2;

      revolver[0].check_filter[count_read] = 0 ;

      revolver[0].count_read ++;
      
      if ( ! ( count_read % reading_chunk ) ) {
	break;
      }
    }

  revolver[0].count_read = count_read;

  const double walltime1(get_wall_time());
  const double cputime1(get_cpu_time());
  
  read_wall_clock += walltime1 - walltime0;
  read_cpu_time += cputime1 - cputime0;
}



void encoding(size_t kmer_size, int num_threads, int chunk_size){
	      // vector<string> & id_vec1, vector<string> & read_vec1, vector<string> & qual_vec1, 
	      // vector<string> & id_vec2, vector<string> & read_vec2, vector<string> & qual_vec2,
	      // vector<bool> & check_filter, size_t & count_read){
  const double walltime0(get_wall_time());
  const double cputime0(get_cpu_time());

  /*
  omp_set_num_threads(num_threads);

#pragma omp parallel for schedule(static,chunk_size) 
  for ( size_t I = 0; I < revolver[1].count_read ; I++)
    revolver[1].check_filter[I] = search_kmer(revolver[1].read_vec1[I], revolver[1].qual_vec1[I], revolver[1].read_vec2[I], revolver[1].qual_vec2[I], kmer_size);
  */

  thread pool[num_threads];

  for ( int i = 0 ; i < num_threads ; i ++ ){
    size_t start = i * chunk_size;
    size_t end = start + chunk_size;
    pool[i] = thread(find_filtered_read,start,end,kmer_size);
  }

  for ( int i = 0 ; i < num_threads ; i ++ ){
    // string out= "Joined: " + std::to_string(i) + "\n";
    // std::cout << out ;
    pool[i].join();
  }

  const double walltime1(get_wall_time());
  const double cputime1(get_cpu_time());
  
  encoding_wall_clock += walltime1 - walltime0;
  encoding_cpu_time += cputime1 - cputime0;
}



void write_fastq(){
// vector<string> & id_vec1, vector<string> & read_vec1, vector<string> & qual_vec1, 
// 		 vector<string> & id_vec2, vector<string> & read_vec2, vector<string> & qual_vec2,
// 		 vector<bool> & check_filter, size_t & count_read){
  const double walltime0(get_wall_time());
  const double cputime0(get_cpu_time());
  

  string output1;
  string output2;
  for ( size_t I = 0 ; I < revolver[2].count_read ; I++)
    {
      output1 = revolver[2].id_vec1[I] + "\n" + revolver[2].read_vec1[I] + "\n+\n" + revolver[2].qual_vec1[I] + "\n" ;
      output2 = revolver[2].id_vec2[I] + "\n" + revolver[2].read_vec2[I] + "\n+\n" + revolver[2].qual_vec2[I] + "\n" ;
      if ( revolver[2].check_filter[I] ) fwrite (output1.c_str(), sizeof(char), output1.size(), outfile_1_out) ;
      if ( revolver[2].check_filter[I] ) fwrite (output2.c_str(), sizeof(char), output2.size(), outfile_2_out) ;
    } // for end

  const double walltime1(get_wall_time());
  const double cputime1(get_cpu_time());
  
  write_wall_clock += walltime1 - walltime0;
  write_cpu_time += cputime1 - cputime0;
}

template < class T >
inline void swap_vec(vector<T> & vec1, vector<T> & vec2){
  T tmp;
  size_t Size = vec1.size();
  for ( size_t i = 0 ; i < Size ; i ++ ) {
    tmp = vec1[i];
    vec1[i] = vec2[i];
    vec2[i] = tmp ;
  }
}

template < class T >
inline void swap_num( T & a, T & b){
  T tmp;
  tmp = a;
  a = b;
  b = tmp ;
}



bool search_kmer(string & read1, string & qual1, size_t kmer_size)
{
  size_t checker1(0);

  uint64_t encode_k(0);
  uint64_t reversed(0);

  uint64_t canonical;
  
  size_t Ncount(0) ;
  size_t i;

  /////////////////////////////////////////////////////////////////////////
  //
  // count non-ref. k-mers in read_1
  //
  
  if ( read1.size() >= kmer_size ){

    for ( i = 0 ; i < kmer_size ; i ++ ){
      encode_forw(encode_k,read1[i],kmer_size);
      encode_reve(reversed,read1[i],kmer_size);
      Ncount += read1[i] == 'N' ;
    }
    
    if ( Ncount == 0 ) {
      canonical = select_canonical( encode_k, reversed ) ;
      if ( kmerDB_set.find(canonical) != kmerDB_set.end() ) checker1 ++ ;
    }
    
    for ( ; i < read1.length() ; i++) {
      encode_forw(encode_k,read1[i],kmer_size);
      encode_reve(reversed,read1[i],kmer_size);
      Ncount += read1 [ i ] == 'N' ;
      Ncount -= read1 [ i - kmer_size ] == 'N' ;
      
      if ( Ncount == 0 ) {
	canonical = select_canonical( encode_k, reversed ) ;
	if ( kmerDB_set.find(canonical) != kmerDB_set.end() ) checker1 ++ ;
      }
    }
  }

  // return result
  return ( checker1 >= non_ref_kmer_cut ) ;
}






void sam_worker(string infile, string outfile, int kmer_size){
  std::ifstream fin ( infile.c_str() );
  std::ofstream fout ( outfile.c_str() );

  vector < string > sam_vec;
  sam_vec.resize(11);
  string tags;

  std::string read;
  std::string qual;

  while ( fin >> sam_vec[0] )
    {
      
      for ( size_t i = 1 ; i < 11 ; i ++ ){
	fin >> sam_vec[i] ;
      }
      read = sam_vec[9];
      qual = sam_vec[10];
      std::getline ( fin , tags );
      if ( search_kmer(read, qual, (unsigned) kmer_size) ){
	fout << sam_vec[0] ;
	for ( size_t i = 1 ; i < 11 ; i ++ ){
	  fout << "\t" << sam_vec[i] ;
	}
	fout << tags << "\n";
      }
    }
  remove ( infile.c_str() );
  fin.close();
  fout.close();
}









int collector_bam (string kmer_table, int kmer_size, string out_prefix, string read_file, int NT)
{
  cout << "[Collector start]\n" ;

  int num_threads = NT;
  std::string tmp;

  const double walltime0(get_wall_time());
  const double cputime0(get_cpu_time());

  // Load k-mer set
  //cout << "Loading k-mers\n";
  float mlf = 0.1 ;

  kmerDB_set.max_load_factor(mlf);
  string kmer;
  string kmer_count;
  
  ifstream fin(kmer_table.c_str());
  while(fin>>kmer>>kmer_count){
    kmerDB_set.insert(return_encode(kmer));
  }
  fin.close();
  
  const double walltime1(get_wall_time());
  const double cputime1(get_cpu_time());
  cout << "[K-mer loading complete] " << walltime1 - walltime0 << " sec (wall-clock), " << cputime1 - cputime0 << " sec (CPU)" << endl;
  
  string outfile_name = out_prefix + ".bam";

  FILE* outfile_bam_out;

  if ( ( outfile_bam_out = fopen (outfile_name.c_str(),"wb" ) ) == NULL ) {
    std::cout << "ERROR!!! " << outfile_name << "cannot be opened properly.\n";
    return EXIT_FAILURE; 
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::string sanm = read_file + "_";
  std::string header = read_file + ".header";
  std::string list_file ;
  list_file = sanm + "list";

  std::string command ;

  command = "samtools view -H " + read_file + " > " + header;
  std::cout << command << "\n" ; 
  system(command.c_str());

  thread sa_tag;
  thread nm_tag;
  std::string command1 = "samtools view -@ " + std::to_string (num_threads/2) + " " + read_file + " | grep \"SA:Z:\" > " + read_file + "_SA";
  std::cout << command1 << "\n" ; 

  std::string command2 =  "samtools view -@ " + std::to_string (num_threads/2) + " " + read_file + " | grep -v \"SA:Z:\" | grep -v \"NM:i:0\" > " + read_file + "_NM";
  std::cout << command2 << "\n" ; 

  sa_tag=thread(system,command1.c_str());
  nm_tag=thread(system,command2.c_str());

  sa_tag.join();
  nm_tag.join();

  command = "cat " + read_file + "_SA " + read_file + "_NM > " + sanm;
  std::cout << command << "\n" ; 
  system(command.c_str());

  tmp = read_file + "_SA" ;
  remove ( tmp.c_str() );
  tmp = read_file + "_NM" ;
  remove ( tmp.c_str() );

  command = "split -n l/" + std::to_string(num_threads) + " --numeric-suffixes=1 " + sanm + " " + sanm;
  std::cout << command << "\n" ; 
  system(command.c_str());

  command = "rm -f " + sanm;
  std::cout << command << "\n" ; 
  system(command.c_str());

  command = "ls " + sanm + "* | grep -v " + read_file + "_SA | grep -v " + read_file + "_NM | grep -v " + list_file + " > " + list_file;
  std::cout << command << "\n" ; 
  system(command.c_str());


  const double walltime2(get_wall_time());
  const double cputime2(get_cpu_time());
  cout << "[Split complete] " << walltime2 - walltime1 << " sec (wall-clock), " << cputime2 - cputime1 << " sec (CPU)" << endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  std::vector < std::string > infile_vec;
  fin.open ( list_file.c_str() );
  while ( fin >> tmp ){
    infile_vec.push_back(tmp);
  }
  fin.close();

  remove ( list_file.c_str() );

  std::vector < std::string > outfile_vec;
  for ( int i = 0 ; i < num_threads ; i ++ ){
    tmp = out_prefix + "_" + std::to_string(i);
    outfile_vec.push_back(tmp);
  }

  //omp_set_num_threads(num_threads);

// #pragma omp parallel
//   {
//     std::string infile ;
//     int tid = omp_get_thread_num();
//     sam_worker(infile_vec[tid], outfile_vec[tid], num_threads);
//   }

  thread pool[num_threads];
  for ( int i = 0 ; i < num_threads ; i ++ ){
    pool[i] = thread ( sam_worker, infile_vec[i], outfile_vec[i], kmer_size);
  }

  for ( int i = 0 ; i < num_threads ; i ++ ){
    pool[i].join();
  }

  command = "cat " + header + " " + out_prefix + "_* | samtools view -@ 30 -Sbh - > " + outfile_name;
  std::cout << command << "\n" ; 
  system ( command.c_str() );
  remove ( header.c_str() );
  for ( int i = 0 ; i < num_threads ; i ++ ){
    remove ( outfile_vec[i].c_str() );
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////


  const double walltime3(get_wall_time());
  const double cputime3(get_cpu_time());

  cout << "[Filtering complete] " << walltime3 - walltime2 << " sec (wall-clock), " << cputime3 - cputime2 << " sec (CPU)" << endl;

  fclose(outfile_bam_out);

  const double walltime4(get_wall_time());
  const double cputime4(get_cpu_time());
  cout << "[read_collector complete]" << endl;
  cout << "Total wall-clock time: " << walltime4 - walltime0 << " sec" << endl;
  cout << "Total CPU time: " << cputime4 - cputime0 << " sec" << endl;
  cout << "Peak memory: " << get_peak_memory() << " KB" << endl; 
  
  return 0;
}


int collector (string kmer_table, int kmer_size, string out_prefix, string read_file_1, string read_file_2, int gz_check, int NT)
{
  cout << "[Collector start]\n" ;

  int num_threads = NT;

  const double walltime0(get_wall_time());
  const double cputime0(get_cpu_time());

  read_wall_clock=0;
  read_cpu_time=0;

  encoding_wall_clock=0;
  encoding_cpu_time=0;

  write_wall_clock=0;
  write_cpu_time=0;

  num_threads -= 2;
  if ( num_threads < 1 ) num_threads = 1;
  
  // Load k-mer set
  //cout << "Loading k-mers\n";
  float mlf = 0.1 ;

  kmerDB_set.max_load_factor(mlf);
  string kmer;
  string kmer_count;
  
  ifstream fin(kmer_table.c_str());
  while(fin>>kmer>>kmer_count){
    kmerDB_set.insert(return_encode(kmer));
  }
  fin.close();
  
  const double walltime1(get_wall_time());
  const double cputime1(get_cpu_time());
  cout << "[K-mer loading complete] " << walltime1 - walltime0 << " sec (wall-clock), " << cputime1 - cputime0 << " sec (CPU)" << endl;
  
  string outfile_name_1_out = out_prefix + "_1.fastq";
  string outfile_name_2_out = out_prefix + "_2.fastq";
    
  if ( ( outfile_1_out = fopen (outfile_name_1_out.c_str(),"wb" ) ) == NULL ) return EXIT_FAILURE; 
  if ( ( outfile_2_out = fopen (outfile_name_2_out.c_str(),"wb" ) ) == NULL ) return EXIT_FAILURE; 

  if ( gz_check == 0 ){
    readFile_1 = fopen(read_file_1.c_str(),"r");
    readFile_2 = fopen(read_file_2.c_str(),"r");
  }
  else{
    readFile_1_gz.open(read_file_1.c_str());
    readFile_2_gz.open(read_file_2.c_str());
  }

  // Check if input files exit. If not, this program is terminated.
  if ( gz_check == 0 ){
    if ( readFile_1 == NULL )
      {
	cout << "ERROR: Input fastq file does not exist" << endl;
	return EXIT_FAILURE ;
      }
    if ( readFile_2 == NULL )
      {
	cout << "ERROR: Input fastq file does not exist" << endl;
	return EXIT_FAILURE ;
      }
  }
  else {
    if ( ! readFile_1_gz )
      {
	cout << "ERROR: Input fastq file does not exist" << endl;
	return EXIT_FAILURE ;      
      }
    if ( ! readFile_2_gz )
      {
	cout << "ERROR: Input fastq file does not exist" << endl;
	return EXIT_FAILURE ;
      }
  }


  // Chunk sizes for openMP
  chunk_size=10000;
  reading_chunk=chunk_size * num_threads;

  for ( size_t i = 0 ; i < 4 ; i ++ ){
    revolver[i].init(reading_chunk);
  }

  size_t total_count_read;
    
  total_count_read = 0;
    
  double wall0;
  double cpu0;
  double wall1;
  double cpu1;
    

  thread READING;
  thread ENCODING;
  thread WRITING;

  /////////////////////////////////////////////////////////////////
  //
  // Initialization-1
  //

  wall0=get_wall_time();
  cpu0=get_cpu_time();

  if ( gz_check == 0 ){
    READING=thread(read_fastq);
  }
  else {
    READING=thread(read_fastq_gz);
  }
  READING.join();

  wall1=get_wall_time();
  cpu1=get_cpu_time();
  total_count_read += revolver[0].count_read;
  cout << "[Number of processed read pairs: " << total_count_read << "] ";
  cout << wall1 - wall0 << " sec (wall-clock), " << cpu1 - cpu0 << " sec (CPU)" ;
  cout << endl;

  //SWAP()
  revolver[3]=revolver[1]; 
  revolver[1]=revolver[0]; 
  revolver[0]=revolver[2]; 
  revolver[2]=revolver[3]; 

  /////////////////////////////////////////////////////////////////
  //
  // Initialization-2
  //

  wall0=get_wall_time();
  cpu0=get_cpu_time();

  if ( gz_check == 0 ){
    if ( ( ! feof(readFile_1) ) && ( ! feof(readFile_1) ) ){
      READING=thread(read_fastq);
    }
  }
  else {
    if ( ( ! readFile_1_gz.eof() ) && ( ! readFile_2_gz.eof() ) ){
      READING=thread(read_fastq_gz);
    }
  }
  ENCODING=thread(encoding,kmer_size,num_threads,chunk_size);

  ENCODING.join();
  READING.join();


  wall1=get_wall_time();
  cpu1=get_cpu_time();
  total_count_read += revolver[0].count_read;
  cout << "[Number of processed read pairs: " << total_count_read << "] ";
  cout << wall1 - wall0 << " sec (wall-clock), " << cpu1 - cpu0 << " sec (CPU)" ;
  cout << endl;

  // SWAP();
  revolver[3]=revolver[1]; 
  revolver[1]=revolver[0]; 
  revolver[0]=revolver[2]; 
  revolver[2]=revolver[3]; 

  /////////////////////////////////////////////////////////////////
  //
  // While loop
  //

  if ( gz_check == 0 ){
    while ( ( ! feof(readFile_1) ) && ( ! feof(readFile_1) ) ){
      
      wall0=get_wall_time();
      cpu0=get_cpu_time();
      
      READING =thread(read_fastq);
      ENCODING=thread(encoding,kmer_size,num_threads,chunk_size);
      WRITING =thread(write_fastq);

		      
      READING.join();
      ENCODING.join();
      WRITING .join();
      
      wall1=get_wall_time();
      cpu1=get_cpu_time();
      total_count_read += revolver[0].count_read;
      cout << "[Number of processed read pairs: " << total_count_read << "] ";
      cout << wall1 - wall0 << " sec (wall-clock), " << cpu1 - cpu0 << " sec (CPU)" ;
      cout << endl;

      // SWAP();      
      revolver[3]=revolver[1]; 
      revolver[1]=revolver[0]; 
      revolver[0]=revolver[2]; 
      revolver[2]=revolver[3]; 
    }
  }
  else{
    while ( ( ! readFile_1_gz.eof() ) && ( ! readFile_1_gz.eof() ) ){
      
      wall0=get_wall_time();
      cpu0=get_cpu_time();
      
      READING=thread(read_fastq_gz);
      ENCODING=thread(encoding,kmer_size,num_threads,chunk_size);
      WRITING =thread(write_fastq);

      READING .join();
      ENCODING.join();
      WRITING .join();
      
      wall1=get_wall_time();
      cpu1=get_cpu_time();

      total_count_read += revolver[0].count_read;

      cout << "[Number of processed read pairs: " << total_count_read << "] ";
      cout << wall1 - wall0 << " sec (wall-clock), " << cpu1 - cpu0 << " sec (CPU)" ;
      cout << endl;

      // SWAP();      
      revolver[3]=revolver[1]; 
      revolver[1]=revolver[0]; 
      revolver[0]=revolver[2]; 
      revolver[2]=revolver[3]; 
    }  
  }

  /////////////////////////////////////////////////////////////////
  //
  // Termination-1
  //
      
  ENCODING=thread(encoding,kmer_size,num_threads,chunk_size);
  WRITING =thread(write_fastq);
  
  ENCODING.join();
  WRITING .join();
    
  wall1=get_wall_time();
  cpu1=get_cpu_time();

  // SWAP();
  revolver[3]=revolver[1]; 
  revolver[1]=revolver[0]; 
  revolver[0]=revolver[2]; 
  revolver[2]=revolver[3]; 

  /////////////////////////////////////////////////////////////////
  //
  // Termination-2
  //

  WRITING =thread(write_fastq);
  WRITING .join();
    
  wall1=get_wall_time();
  cpu1=get_cpu_time();


  /////////////////////////////////////////////////////////////////
  //
  // Finishing
  //

  const double walltime2(get_wall_time());
  const double cputime2(get_cpu_time());

  cout << "\n";
  cout << "[Filtering complete] " << walltime2 - walltime1 << " sec (wall-clock), " << cputime2 - cputime1 << " sec (CPU)" << endl;
  cout << "\tReading:    " << read_wall_clock << " sec (wall-clock), " << read_cpu_time << " sec (cpu time)" << endl;
  cout << "\tCollecting: " << encoding_wall_clock << " sec (wall-clock), " << encoding_cpu_time << " sec (cpu time)" << endl;
  cout << "\tWriting:    " << write_wall_clock << " sec (wall-clock), " << write_cpu_time << " sec (cpu time)" << endl;
  cout << "\n";
  
  if ( gz_check == 0 ){
    fclose(readFile_1);
    fclose(readFile_2);
  }
  else{
    readFile_1_gz.close();
    readFile_2_gz.close();
  }

  fclose(outfile_1_out);
  fclose(outfile_2_out);

  const double walltime3(get_wall_time());
  const double cputime3(get_cpu_time());
  cout << endl << "[read_collector complete]" << endl;
  cout << "Total wall-clock time: " << walltime3 - walltime0 << " sec" << endl;
  cout << "Total CPU time: " << cputime3 - cputime0 << " sec" << endl;
  cout << "Peak memory: " << get_peak_memory() << " KB" << endl; 

  return 0;
}
