#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <functional>



uint32_t adler32(std::string & input) 
{
  //uint32_t MOD_ADLER = 65521;
  uint32_t a = 1;
  uint32_t b = 0;

  std::size_t Size = input.size();

  for (std::size_t i = 0; i < Size ; i++ ){
    a = (a + input[i]) % 65521;
    b = (b + a) % 65521;
  }
    
  return (b << 16) | a;
}


int64_t file_size(std::string infile){
  std::ifstream fin;
  fin.open(infile.c_str(), std::ios_base::binary);
  fin.seekg(0,std::ios_base::end);
  int64_t size = fin.tellg();
  fin.close();
  return size;
}



void sort_and_write_fastq ( std::string infile , std::string outfile ){

  std::ifstream fin ( infile.c_str() );
  std::ofstream fout ( outfile.c_str() );

  std::string id;
  std::string seq;
  std::string desc;
  std::string qual;

  std::map < std::string , std::pair < std::string , std::string > > fastq_map;

  while ( std::getline ( fin , id ) && std::getline ( fin , seq ) && 
	  std::getline ( fin , desc ) && std::getline ( fin , qual ) ){
    if ( id.substr(id.size()-2) == "\1" || id.substr(id.size()-2) == "\2"){
      id = id.substr(0,id.size()-2);
    }
    fastq_map[id].first = seq;
    fastq_map[id].second = qual;
  }

  for ( auto & i : fastq_map ){
    id = i.first;
    seq = i.second.first;
    qual = i.second.second;
    fout << id << "\n" << seq << "\n+\n" << qual << "\n";
  }

  fin.close();
  fout.close();
}




void split_file ( std::string infile , std::vector < std::string > fname_vec, int64_t fsize ){
  std::hash<std::string> hashing;
  std::ifstream fin;
  std::string id;
  std::string seq;
  std::string desc;
  std::string qual;

  // open file streams
  fin.open(infile.c_str());
  std::vector < std::ofstream > fout;
  fout.resize(fsize);
  for ( uint i = 0 ; i < fsize ; i ++ ){
    fout[i].open(fname_vec[i].c_str());
  }


  while ( std::getline ( fin , id ) && std::getline ( fin , seq ) && 
	  std::getline ( fin , desc ) && std::getline ( fin , qual ) ){
    if ( id.substr(id.size()-2) == "\1" || id.substr(id.size()-2) == "\2"){
      id = id.substr(0,id.size()-2);
    }
    std::size_t hval = hashing ( id );
    // std::size_t hval = adler32 ( id );
    std::size_t fnumber = hval % fsize;
    fout[fnumber] << id << "\n" << seq << "\n+\n" << qual << "\n";
  }

  // close file streams
  fin.close();
  for ( uint i = 0 ; i < fsize ; i ++ ){
    fout[i].close();
  }

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Main
//

int main ( int argc , char ** argv ){
  if ( argc == 1 ){
    std::cout << "sort_fastq input.fastq output.fastq chunk_file_size(GB)\n";
    return 0;
  }

  std::string infile(argv[1]);
  std::string outfile(argv[2]);
  int64_t FSIZE = atoi(argv[3]);

  int64_t GB = 1024*1024*1024 ;
  int64_t fsize = file_size(infile) / ( FSIZE * GB ) + 1;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // If file size is < 1 GB, just sort using binary tree (std::map).
  //

  //std::cout << "File size: " << file_size(infile) << " Byte\n";

  if ( fsize == 1 ){
    //std::cout << "Just sort\n";
    sort_and_write_fastq(infile, outfile);
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // If file size is greater 1 GB, split files 
  // 
  // < 2 GB --> 2 files
  // < 3 GB --> 3 files
  // .
  // < n GB --> n files
  else{
    // std::cout << "Split sort\n";

    std::vector < std::string > fname_vec;
    std::vector < std::string > sname_vec;

    // initializing variables
    for ( uint i = 0 ; i < fsize ; i ++ ){
      std::string fname = outfile + "_" + std::to_string(i);
      std::string sname = fname + "_sort";
      sname_vec.push_back(sname);
      fname_vec.push_back(fname);
    }

    // split file
    split_file ( infile, fname_vec, fsize );

    // sort by id
    for ( uint i = 0 ; i < fsize ; i ++ ){
      std::string fname = fname_vec[i];
      std::string sname = sname_vec[i];
      sort_and_write_fastq(fname, sname);
    }

    // remove split files
    for ( uint i = 0 ; i < fsize ; i ++ ){
      remove(fname_vec[i].c_str());
    }

    // merge and print64_t
    std::ofstream fout ( outfile.c_str() );
    for ( uint i = 0 ; i < fsize ; i ++ ){
      std::string sname = sname_vec[i];
      std::string id;
      std::string seq;
      std::string desc;
      std::string qual;
      std::ifstream fin ( sname.c_str() );
      while ( std::getline ( fin , id ) && std::getline ( fin , seq ) &&
	      std::getline ( fin , desc ) && std::getline ( fin , qual ) ){
	// if ( id.substr(id.size()-2) == "\1" || id.substr(id.size()-2) == "\2"){
	//   id = id.substr(0,id.size()-2);
	// }
	fout << id << "\n" << seq << "\n+\n" << qual << "\n";
      }
      fin.close();
    }

    // remove sorted split files
    for ( uint i = 0 ; i < fsize ; i ++ ){
      remove(sname_vec[i].c_str());
    }

    fout.close();
  }

  return 0;
}
