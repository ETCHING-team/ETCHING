//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------


#include "functions.hpp" 



int check_concord ( const Position & Pos1, const Position & Pos2, const int & insert_size ){
  return Pos1.first == Pos2.first && ( Pos1.second - Pos2.second < 2 * insert_size && Pos2.second - Pos1.second < 2 * insert_size) ? 1 : 0 ;
}



int FindConcordant ( const Position & Pos_read , const Position & Pos_split , const Position & Pos_mate , const int & insert_size ){
  return check_concord(Pos_read,Pos_split,insert_size) + check_concord(Pos_read,Pos_mate,insert_size) + check_concord(Pos_mate,Pos_split,insert_size) == 1 ? 1 : 0 ;
}



std::vector < BamTools::CigarOp > cigar_parse ( std::string cigar ){
  std::vector < BamTools::CigarOp > output;

  std::string Type;
  int Length;
  BamTools::CigarOp op;

  while ( cigar.size() > 0 ){
    Length = atoi(cigar.c_str());
    std::size_t found = std::to_string(Length).size();
    Type=cigar.substr(found,1).c_str();
    op.Type = Type[0];
    op.Length = Length;
    output.push_back(op);
    cigar = cigar.substr(found+1);
  }

  return output;
}




std::string parse_SA_tag ( std::string & Tag ){
  std::string output;
  std::size_t found = Tag.find(";");
  output = Tag.substr(0,found);
  Tag = Tag.substr(found +1 );
  return output;
}



Position find_split ( const std::string & tag , const std::unordered_map < std::string , int > & id_ref_map ){
  int id_num;
  int pos;

  std::size_t found;
  std::string pref;
  std::string suff;
  std::string cigar;

  std::vector < BamTools::CigarOp > cigar_vector;

  // finding chr id
  found = tag.find(",");
  pref = tag.substr(0,found);
  id_num = id_ref_map.find(pref)->second;

  // finding position
  suff = tag.substr(found+1);
  found = suff.find(",");
  pref = suff.substr(0,found);
  pos = atoi ( pref.c_str() ) - 1;

  // finding cigar
  suff = suff.substr(found+3);
  found = suff.find(",");
  cigar = suff.substr(0,found);
  cigar_vector = cigar_parse(cigar);

  int size = cigar_vector.size();

  bool rev_condition = ( cigar_vector[0].Type == 'S' || cigar_vector[0].Type == 'H' ) ;
  bool for_condition = ( cigar_vector[size-1].Type == 'S' || cigar_vector[size-1].Type == 'H' ) ;

  //                      /
  //                     /
  //                    /
  //   ======D=========/
  //       delet       Clpped --> calculating shift
  //

  char Type;
  int Length;

  if ( rev_condition == 0 ){
    if (for_condition == 1 ){
      for ( int i = 0 ; i < size - 1 ; i ++ ){
	Type = cigar_vector[i].Type;
	Length = cigar_vector[i].Length;
	if ( Type == 'D' || Type == 'M' || Type == '=' || Type == 'X' ){
	  pos += Length ;
	}
      }
      pos --;
    }
    else {
      pos = -1;
    }
  }
  else if ( for_condition == 1 ){
    pos = - 1 ;
  }

  return std::make_pair (id_num,pos);
}



std::string mate_strand ( const std::string & tag ) {

  std::size_t found;
  std::string pref;
  std::string suff;
  std::string cigar;

  std::vector < BamTools::CigarOp > cigar_vector;

  // finding chr id
  found = tag.find(",");
  pref = tag.substr(0,found);

  // finding position
  suff = tag.substr(found+1);
  found = suff.find(",");
  pref = suff.substr(0,found);

  // finding cigar
  suff = suff.substr(found+3);
  found = suff.find(",");
  cigar = suff.substr(0,found);
  cigar_vector = cigar_parse(cigar);

  int size = cigar_vector.size();

  if ( ( cigar_vector[0].Type == 'S' || cigar_vector[0].Type == 'H' ) &&
       !( cigar_vector[size-1].Type == 'S' || cigar_vector[size-1].Type == 'H' ) ) {
    return "R";
  }

  if ( !( cigar_vector[0].Type == 'S' || cigar_vector[0].Type == 'H' ) &&
       ( cigar_vector[size-1].Type == 'S' || cigar_vector[size-1].Type == 'H' ) ) {
    return "F";
  }

  return "";
}




std::string currentDate() {
  time_t     now = time(0);
  struct tm  tstruct;
  char       buf[80];
  tstruct = *localtime(&now);
  strftime(buf, sizeof(buf), "%Y-%m-%d", &tstruct); // YYYY-MM-DD

  return buf;
}


std::string currentDateTime() {
  time_t     now = time(0);
  struct tm  tstruct;
  char       buf[80];
  tstruct = *localtime(&now);
  strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct); // YYYY-MM-DD.HH:mm:ss

  return buf;
}
