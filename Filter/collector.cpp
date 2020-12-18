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
  cout << "[K-mer loading complete]" << walltime1 - walltime0 << " sec (wall-clock), " << cputime1 - cputime0 << " sec (CPU)" << endl;
  
  string outfile_name_1_out = out_prefix + "_1.fastq";
  string outfile_name_2_out = out_prefix + "_2.fastq";
    
  if ( ( outfile_1_out = fopen (outfile_name_1_out.c_str(),"wb" ) ) == NULL ) return EXIT_FAILURE; 
  if ( ( outfile_2_out = fopen (outfile_name_2_out.c_str(),"wb" ) ) == NULL ) return EXIT_FAILURE; 

  // Include the kmer list to set container
  // string id1, read1, desc1, qual1;
  // string id2, read2, desc2, qual2;

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
  chunk_size=100;
  fastq_chunk=10000 * num_threads;

  for ( size_t i = 0 ; i < 4 ; i ++ ){

    id_vec1[i].reserve(fastq_chunk);
    read_vec1[i].reserve(fastq_chunk);
    qual_vec1[i].reserve(fastq_chunk);

    id_vec2[i].reserve(fastq_chunk);
    read_vec2[i].reserve(fastq_chunk);
    qual_vec2[i].reserve(fastq_chunk);

    check_filter[i].reserve(fastq_chunk); 
  }
    
  for ( size_t i = 0 ; i < 4 ; i ++ ){

    id_vec1[i].resize(fastq_chunk);
    read_vec1[i].resize(fastq_chunk);
    qual_vec1[i].resize(fastq_chunk);

    id_vec2[i].resize(fastq_chunk);
    read_vec2[i].resize(fastq_chunk);
    qual_vec2[i].resize(fastq_chunk);

    check_filter[i].resize(fastq_chunk); 

    count_read[i]=0;
  }
    
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
  total_count_read += count_read[0];
  cout << "[Number of processed read pairs: " << total_count_read << "] ";
  cout << wall1 - wall0 << " sec (wall-clock), " << cpu1 - cpu0 << " sec (CPU)" ;
  cout << endl;

  SWAP();

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
  ENCODING=thread(encoding,kmer_size,num_threads);
  ENCODING.join();
  READING.join();


  wall1=get_wall_time();
  cpu1=get_cpu_time();
  total_count_read += count_read[0];
  cout << "[Number of processed read pairs: " << total_count_read << "] ";
  cout << wall1 - wall0 << " sec (wall-clock), " << cpu1 - cpu0 << " sec (CPU)" ;
  cout << endl;

  SWAP();
  /////////////////////////////////////////////////////////////////
  //
  // While loop
  //

  if ( gz_check == 0 ){
    while ( ( ! feof(readFile_1) ) && ( ! feof(readFile_1) ) ){
      
      wall0=get_wall_time();
      cpu0=get_cpu_time();
      
      READING=thread(read_fastq);
      ENCODING=thread(encoding,kmer_size,num_threads);
      WRITING =thread(write_fastq);

      ENCODING.join();
      READING .join();
      WRITING .join();
      
      wall1=get_wall_time();
      cpu1=get_cpu_time();
      total_count_read += count_read[0];
      cout << "[Number of processed read pairs: " << total_count_read << "] ";
      cout << wall1 - wall0 << " sec (wall-clock), " << cpu1 - cpu0 << " sec (CPU)" ;
      cout << endl;

      SWAP();      
    }
  }
  else{
    while ( ( ! readFile_1_gz.eof() ) && ( ! readFile_1_gz.eof() ) ){
      
      wall0=get_wall_time();
      cpu0=get_cpu_time();
      
      READING=thread(read_fastq_gz);
      ENCODING=thread(encoding,kmer_size,num_threads);
      WRITING =thread(write_fastq);

      ENCODING.join();
      READING .join();
      WRITING .join();
      
      wall1=get_wall_time();
      cpu1=get_cpu_time();

      total_count_read += count_read[0];

      cout << "[Number of processed read pairs: " << total_count_read << "] ";
      cout << wall1 - wall0 << " sec (wall-clock), " << cpu1 - cpu0 << " sec (CPU)" ;
      cout << endl;

      SWAP();      
    }  
  }

  /////////////////////////////////////////////////////////////////
  //
  // Termination-1
  //
      
  ENCODING=thread(encoding,kmer_size,num_threads);
  WRITING =thread(write_fastq);

  ENCODING.join();
  WRITING .join();
    
  wall1=get_wall_time();
  cpu1=get_cpu_time();

  SWAP();
    
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


  for ( size_t i = 0 ; i < 4 ; i ++ ){

    id_vec1[i].clear();
    read_vec1[i].clear();
    qual_vec1[i].clear();

    id_vec2[i].clear();
    read_vec2[i].clear();
    qual_vec2[i].clear();

    al_vec[i].clear();

    check_filter[i].clear();

  }

  
  return 0;
}
