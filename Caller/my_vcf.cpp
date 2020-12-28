//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------


#include "my_vcf.hpp" 




VCF::VCF(std::string line){
  parse(line);
}



void VCF::parse_etching(std::string line){
  // std::cout << line << "\n";
  line = parse(line);
  // std::cout << line << "\n\n";

  /////////////////////////////////////////////////////////////////
  //
  // Reading FORMAT information
  //
  // feature_str = "CR:SR:PE:MQ:DEPDIF:NXA:TCB:ENTROPY:PURITY:SEQDEP";
  //

  std::size_t found;
  found = line.find('\t');
  info = line.substr(0,found);
  line = line.substr(found+1);
  // std::cout << line << "\n";
  // std::cout << info << "\n\n";

  found = line.find('\t');
  feature_str  = line.substr(0,found);
  line = line.substr(found+1);
  // std::cout << line << "\n";
  // std::cout << feature_str << "\n\n";

  found = line.find('\t');
  feature_str1 = line.substr(0,found);
  line = line.substr(found+1);
  // std::cout << line << "\n";
  // std::cout << feature_str1 << "\n\n";

  found = line.find('\t');
  feature_str2 = line.substr(0,found);
  line = line.substr(found+1);
  // std::cout << line << "\n";
  // std::cout << feature_str2 << "\n\n";

  feature = parse_feature(feature_str1);
  feature1 = parse_feature(feature_str1);
  feature2 = parse_feature(feature_str2);

  input_features ( feature_str1 );
  input_features_2 ( feature_str2 );
 
}


std::string VCF::parse(std::string line){

  std::string tmp;
  std::string key;
  std::size_t sz;
  std::size_t found;
  
  chr2 = "";
  pos2 = -1;


  // chr
  found = line.find('\t');
  chr1 = line.substr(0,found);
  line = line.substr(found+1);
  
  // pos
  found = line.find('\t');
  pos1 =  atoi(line.substr(0,found).c_str());
  line = line.substr(found+1);

  // id
  found = line.find('\t');
  sv_id =  line.substr(0,found);
  line = line.substr(found+1);

  // reference sequence
  found = line.find('\t');
  ref =  line.substr(0,found);
  line = line.substr(found+1);

  // alt sequence
  found = line.find('\t');
  alt =  line.substr(0,found);
  line = line.substr(found+1);

  // quality
  found = line.find('\t');
  if ( isdigit(line.substr(0,found)[0] ) ){
    qual =  atoi(line.substr(0,found).c_str());
  }
  else{
    qual = -1;
  }
  line = line.substr(found+1);

  // filter
  found = line.find('\t');
  filter =  line.substr(0,found);
  line = line.substr(found+1);

  // information
  info =  line;

  /////////////////////////////////////////////////////////////////
  //
  // read chr2 and pos2
  //

  std::string chr2_key="CHR2=";
  std::string pos2_key=";END=";

  tmp = alt;
  std::string tmp1;

  key = "[";
  tmp1 = cut_str(tmp, key);
  if ( tmp.size() > 0 ){
    chr2 = cut_str(tmp,":");
    tmp = cut_str(tmp,"[");
    pos2 = std::stol(tmp);
  }

  tmp = alt;
  key = "]";
  tmp1 = cut_str(tmp, key);
  if ( tmp.size() > 0 ){ 
    chr2 = cut_str(tmp, ":");
    tmp = cut_str(tmp,key);
    pos2 = std::stol(tmp);
  }


  if ( chr2 == "" ){
    tmp = info;
    tmp1 = cut_str(tmp,chr2_key);
    if ( tmp.size() > 0 ){
      key = ";";
      tmp1 = cut_str(tmp,key);
      chr2 = tmp1;
    }
    else{
      chr2 = chr1;
    }
  }

  if ( pos2 == -1 ){
    tmp = info ;
    tmp1 = cut_str(tmp,pos2_key);
    if ( tmp.size() > 0 ){
      pos2 = std::stol(tmp,&sz);
    }
  }

  if ( pos2 == -1 ){
    tmp = info ;
    tmp1 = cut_str(tmp,"END=");
    if ( tmp.size() > 0 ){
      pos2 = std::stol(tmp,&sz);
    }
  }

  /////////////////////////////////////////////////////////////////
  //
  // Reading SVTYPE
  //
  
  tmp = info;
  key = "SVTYPE=";
  cut_str(tmp,key);
  svtype = tmp.substr(0,3);

  if ( svtype == "TRA" ) svtype = "BND";
  

  /////////////////////////////////////////////////////////////////
  //
  // Reading STRAND
  //

  strand = "";
  
  key = "SVTYPE=BND";
  if ( info.find(key) != std::string::npos ){
    if ( alt[alt.size()-1] == '[' ) strand = "FR"; // fixed
    else if ( alt[alt.size()-1] == ']' ) strand = "FF"; // fixed
    else if ( alt[0] == '[' ) strand = "RR"; // fixed
    else if ( alt[0] == ']' ) strand ="RF"; // fixed
  }
  
  if ( alt[alt.size()-1] == '.' ){
    chr2="";
    pos2=-1;
    strand ="F";
  }
  else if ( alt[0] == '.' ){
    chr2="";
    pos2=-1;
    strand ="R";
  }

  if (strand.size() == 0){
    key="STRANDS=";
    tmp=info;
    cut_str ( tmp , key );
    if ( tmp.size() != 0 ){
      tmp = tmp.substr(0,2);
      if ( tmp == "+-" ) strand = "FR"; // fixed
      else if ( tmp == "++" ) strand = "FF"; // fixed
      else if ( tmp == "--" ) strand = "RR"; // fixed
      else if ( tmp == "-+" ) strand = "RF"; // fixed
    }
    else{
      key="CT=";
      tmp=info;
      cut_str(tmp,key);
      tmp=tmp.substr(0,4);
      if ( tmp == "5to3" ){
	strand = "FR"; // fixed
      }
      else if ( tmp == "5to5" ){ 
	strand = "FF"; // fixed
      }
      else if ( tmp == "3to3" ){
	strand = "RR"; // fixed
      }
      else if ( tmp == "3to5" ){
	strand = "RF"; // fixed
      }
    }
  }



  /////////////////////////////////////////////////////////////////
  //
  // MATEID
  //

  tmp = info ;
  key = "MATEID=";
  cut_str(tmp,key);
  mate_id = cut_str(tmp,";");

  return line;
}




Feature VCF::parse_feature(std::string feature_str){
  Feature _feature;

  std::string tmp;
  std::size_t sz;
  std::string key=":";

  tmp = cut_str(feature_str,key);
  double CR = std::stod(tmp,&sz);

  tmp = cut_str(feature_str,key);
  double SR = std::stod(tmp,&sz);
  
  tmp = cut_str(feature_str,key);
  double PE = std::stod(tmp,&sz);
  
  tmp = cut_str(feature_str,key);
  double MQ = std::stod(tmp,&sz);
  
  tmp = cut_str(feature_str,key);
  double DEPDIF = std::stod(tmp,&sz);
  
  tmp = cut_str(feature_str,key);
  double NXA = std::stod(tmp,&sz);
  
  tmp = cut_str(feature_str,key);
  double TCB = std::stod(tmp,&sz);
  
  tmp = cut_str(feature_str,key);
  double ENTROPY = std::stod(tmp,&sz);

  tmp = cut_str(feature_str,key);
  double PURITY = std::stod(tmp,&sz);

  tmp = cut_str(feature_str,key);
  double SEQDEP = std::stod(tmp,&sz);

  
  _feature["CR"]=CR;
  _feature["SR"]=SR;
  _feature["PE"]=PE;
  _feature["MQ"]=MQ;
  _feature["DEPDIF"]=DEPDIF;
  _feature["NXA"]=NXA;
  _feature["TCB"]=TCB;
  _feature["ENTROPY"]=ENTROPY;
  _feature["PURITY"]=PURITY;
  _feature["SEQDEP"]=SEQDEP;

  return _feature;
}


void VCF::input_features ( std::string feature_str ){
  std::string tmp;
  std::size_t sz;
  std::string key=":";

  tmp = cut_str(feature_str,key);
  cr = std::stod(tmp,&sz);

  tmp = cut_str(feature_str,key);
  sr = std::stod(tmp,&sz);

  tmp = cut_str(feature_str,key);
  pe = std::stod(tmp,&sz);

  tmp = cut_str(feature_str,key);
  mq = std::stod(tmp,&sz);

  tmp = cut_str(feature_str,key);
  depdif = std::stod(tmp,&sz);

  tmp = cut_str(feature_str,key);
  nxa = std::stod(tmp,&sz);

  tmp = cut_str(feature_str,key);
  tcb = std::stod(tmp,&sz);

  tmp = cut_str(feature_str,key);
  entropy = std::stod(tmp,&sz);

  tmp = cut_str(feature_str,key);
  purity = std::stod(tmp,&sz);

  tmp = cut_str(feature_str,key);
  seqdep = std::stod(tmp,&sz);
}



void VCF::input_features_2 ( std::string feature_str ){
  std::string tmp;
  std::size_t sz;
  std::string key=":";

  tmp = cut_str(feature_str,key);
  cr2 = std::stod(tmp,&sz);

  tmp = cut_str(feature_str,key);
  sr2 = std::stod(tmp,&sz);

  tmp = cut_str(feature_str,key);
  pe2 = std::stod(tmp,&sz);

  tmp = cut_str(feature_str,key);
  mq2 = std::stod(tmp,&sz);

  tmp = cut_str(feature_str,key);
  depdif2 = std::stod(tmp,&sz);

  tmp = cut_str(feature_str,key);
  nxa2 = std::stod(tmp,&sz);

  tmp = cut_str(feature_str,key);
  tcb2 = std::stod(tmp,&sz);

  tmp = cut_str(feature_str,key);
  entropy2 = std::stod(tmp,&sz);

  tmp = cut_str(feature_str,key);
  purity = std::stod(tmp,&sz);

  tmp = cut_str(feature_str,key);
  seqdep = std::stod(tmp,&sz);
}


std::string VCF::cut_str(std::string & input, std::string key){
  std::size_t found = input.find(key);
  std::string output;
  if ( found != std::string::npos ){
    output = input.substr(0,found);
    input = input.substr(found+key.size());
  }
  else{
    output = input;
    input.clear();
  }
  return output;
}

  



void VCF::modify_svtype_info(){
  std::size_t start,end;
  std::string svt="SVTYPE=";
  std::string front;
  std::string back;
  if ( svtype=="BND"){
    start = info.find(svt) + svt.size();
    end = start + 3;
    front = info.substr(0,start);
    back = info.substr(end);
    info = front + "BND" + back;
  }
  else if ( svtype=="INS"){
    start = info.find(svt) + svt.size();
    end = start + 3;
    front = info.substr(0,start);
    back = info.substr(end);
    info = front + "INS" + back;
  }
  else if ( svtype=="DEL" ){
    start = info.find(svt) + svt.size();
    end = start + 3;
    front = info.substr(0,start);
    back = info.substr(end);
    info = front + "DEL" + back;
  }
  else if ( svtype=="DUP" ){
    start = info.find(svt) + svt.size();
    end = start + 3;
    front = info.substr(0,start);
    back = info.substr(end);
    info = front + "DUP" + back;
  }
  else if ( svtype=="INV" ){
    start = info.find(svt) + svt.size();
    end = start + 3;
    front = info.substr(0,start);
    back = info.substr(end);
    info = front + "INV" + back;
  }
}


void VCF::make_info(){
  info.clear();
  // info += "SVMETHOD=ETCHING;";
  if ( svtype.size() != 0 ) info += "SVTYPE=" + svtype + ";";
  if ( chr2.size() != 0 ) info += "CHR2=" + chr2 + ";";
  if ( pos2 > 0 ) info += "END=" + std::to_string(pos2) + ";";
  if ( strand.size() !=0 ) info += "STRAND=" + strand + ";";
  if ( mate_id.size() !=0 ) {
    info += "MATEID=" + mate_id + ";";
  }
  

  if ( strand == "F" ){
    info +="REPATH="   + chr1 + ":" + std::to_string(pos1) + "(+)-UNKNOWN" + ";";
  }
  else if ( strand == "R" ){
    info +="REPATH="   + chr1 + ":" + std::to_string(pos1) + "(-)-UNKNOWN" + ";";
  }
  else if ( strand == "FR" ){
    info +="REPATH="   + chr1 + ":" + std::to_string(pos1) + "(+)-" + chr2 + ":" + std::to_string(pos2) + "(+)" + ";";
  }
  else if ( strand == "FF" ){
    info +="REPATH="   + chr1 + ":" + std::to_string(pos1) + "(+)-" + chr2 + ":" + std::to_string(pos2) + "(-)" + ";";
  }
  else if ( strand == "RR" ){
    info +="REPATH="   + chr1 + ":" + std::to_string(pos1) + "(-)-" + chr2 + ":" + std::to_string(pos2) + "(+)" + ";";
  }
  else if ( strand == "RF" ){
    info +="REPATH="   + chr1 + ":" + std::to_string(pos1) + "(-)-" + chr2 + ":" + std::to_string(pos2) + "(-)" + ";";
  }


  if ( chr1 == chr2 ){
    svlen = pos1 - pos2;
    if (svlen<0) svlen = -svlen;
    info += "SVLEN=" + std::to_string(svlen) + ";";
  }
  if ( chr1 != chr2 && svtype != "SND" ) info += "INTERCHR;";
  // info += ";SOMATIC";

  if ( info[info.size()-1] == ';' ) info.pop_back();

  feature_str = "CR:SR:PE:MQ:DEPDIF:NXA:TCB:ENTROPY:PURITY:SEQDEP";
  
  
  feature_str1 = std::to_string(cr) ;
  feature_str1 += ":" + std::to_string(sr) ;
  feature_str1 += ":" + std::to_string(pe) ;
  feature_str1 += ":" + std::to_string(mq) ; // Mapping quality
  feature_str1 += ":" + std::to_string(depdif) ;
  feature_str1 += ":" + std::to_string(nxa) ;
  feature_str1 += ":" + std::to_string(tcb) ;
  feature_str1 += ":" + std::to_string(entropy) ;
  feature_str1 += ":" + std::to_string(purity) ; // Tumor purity
  feature_str1 += ":" + std::to_string(seqdep) ; // Sequencing depth

  feature_str2 = std::to_string(cr2) ;
  feature_str2 += ":" + std::to_string(sr2) ;
  feature_str2 += ":" + std::to_string(pe2) ;
  feature_str2 += ":" + std::to_string(mq2) ; // Mapping quality
  feature_str2 += ":" + std::to_string(depdif2) ;
  feature_str2 += ":" + std::to_string(nxa2) ;
  feature_str2 += ":" + std::to_string(tcb2) ;
  feature_str2 += ":" + std::to_string(entropy2) ;
  feature_str2 += ":" + std::to_string(purity) ; // Tumor purity
  feature_str2 += ":" + std::to_string(seqdep) ; // Sequencing depth

}




std::string VCF::to_string(){
  std::string tmp;
  std::string Qual=".";
  if (qual>=0 ) Qual=std::to_string(qual);
  tmp = chr1 + "\t" + std::to_string(pos1) + "\t" + sv_id + "\t" + ref + "\t" + alt + "\t" + Qual + "\t" + filter + "\t" + info ;
  if ( feature_str.size() > 0 ) tmp += "\t" + feature_str;
  if ( feature_str1.size() > 0 ) tmp += "\t" + feature_str1;
  if ( feature_str2.size() > 0 ) tmp += "\t" + feature_str2;
  return tmp;
}


std::string VCF::to_string_short(){
  std::string tmp;
  std::string Qual=".";
  if (qual>0 ) Qual=std::to_string(qual);
  tmp = chr1 + "\t" + std::to_string(pos1) + "\t" + sv_id + "\t" + ref + "\t" + alt + "\t" + Qual + "\t" + filter + "\t" + info ;
  // if ( feature_str.size() > 0 ) tmp += "\t" + feature_str;
  // if ( feature_str1.size() > 0 ) tmp += "\t" + feature_str1;
  // if ( feature_str2.size() > 0 ) tmp += "\t" + feature_str2;
  return tmp;
}


void VCF::resize_tool_comp(int Size){
  tool_comp.clear();
  tool_comp.resize(Size);
}


VCF return_mate(VCF input){
  VCF output = input;

  output.chr1 = input.chr2;
  output.pos1 = input.pos2;

  output.chr2 = input.chr1;
  output.pos2 = input.pos1;

  if ( input.strand == "FR" ) output.strand == "RF";
  else if ( input.strand == "RF" ) output.strand == "FR";

  return output;
}




///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////


void VCF_CLASS::make_header(){
  make_header_short();
  metainfo = metainfo + "\n"
    + "##FILTER=<ID=PASS,Description=\"High confidential variation greater than or equal to cut-off.\">\n"
    + "##FILTER=<ID=LOWQUAL,Description=\"Low quality variation lower than cut-off.\">\n"
    + "##FORMAT=<ID=CR,Number=1,Type=Integer,Description=\"Number of clipped reads supporting the variant.\">\n"
    + "##FORMAT=<ID=SR,Number=1,Type=Integer,Description=\"Number of split reads supporting the variant.\">\n"
    + "##FORMAT=<ID=PE,Number=1,Type=Integer,Description=\"Number of paired-end reads supporting the variant.\">\n"
    + "##FORMAT=<ID=MQ,Number=1,Type=Integer,Description=\"Average mapping quality of the reads supporting the variant.\">\n"
    + "##FORMAT=<ID=DEPDIF,Number=1,Type=Float,Description=\"Depth difference between right and left side of the variant.\">\n"
    + "##FORMAT=<ID=NXA,Number=1,Type=Integer,Description=\"Total number of alternative alignments (XA tag) of the reads supporting the variant.\">\n"
    + "##FORMAT=<ID=TCB,Number=1,Type=Integer,Description=\"Total clipped length of the reads supporting the variant.\">\n"
    + "##FORMAT=<ID=ENTROPY,Number=1,Type=Float,Description=\"Entropy of BP path. Zero means only one path. The more paths, the more entropy.\">\n"
    + "##FORMAT=<ID=PURITY,Number=1,Type=Float,Description=\"Tumor purity of the sample.\">\n"
    + "##FORMAT=<ID=SEQDEP,Number=1,Type=Float,Description=\"Sequencing depth of the sample.\">";
  header = header +"\tFORMAT\tFIRST_MATE\tSECOND_MATE";
}


void VCF_CLASS::make_header_short(){
  metainfo = "##fileformat=VCFv4.2\n";
  metainfo = metainfo
    + "##fileDate="+currentDate()+"\n"
    + "##source=" + etching_version + "\n"
    + "##reference_genome="+reference+"\n";
  for ( auto & i : genome_info ){
    metainfo +="##contig=<ID="+i.first+",length="+std::to_string(i.second)+">\n";
  }
  metainfo = metainfo
    + "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant.\">\n"
    + "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation.\">\n"
    + "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant.\">\n"
    + "##INFO=<ID=STRAND,Number=1,Type=String,Description=\"Strand(s) at break-points, F (5' or +) or R (3' or -).\">\n"
    + "##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate BND, or original IDs of BNDs before SV typing.\">\n"
    + "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV.\">\n"
    // + "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic variation.\">\n"
    // + "##INFO=<ID=GERMLINE,Number=0,Type=Flag,Description=\"Germline variant.\">\n"
    + "##INFO=<ID=REPATH,Number=1,Type=String,Description=\"Path of rearrangement.\">\n"
    + "##INFO=<ID=INTERCHR,Number=0,Type=Flag,Description=\"Inter-chromosomal variation.\">\n"
    // + "##INFO=<ID=SEQTYPE,Number=1,Type=String,Description=\"Sequencing data type. One of WGS, WES and PANEL.\">\n"
    // + "##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Type of approach used to detect SV.\">\n"
    // + "##INFO=<ID=FUSION,Number=.,Type=String,Description=\"Information of fusion gene that is related with the variation.\">\n"
    // + "##ALT=<ID=INS,Description=\"Insertion\">\n"
    + "##ALT=<ID=DEL,Description=\"Deletion\">\n"
    + "##ALT=<ID=DUP,Description=\"Duplication\">\n"
    + "##ALT=<ID=INV,Description=\"Inversion\">\n"
    + "##ALT=<ID=BND,Description=\"Breakend, which may be translocation or unclssified variation\">\n"
    + "##ALT=<ID=SND,Description=\"Single breakend, which may be translocation to unknown contig/scaffold or unclssified variation\">";

  header
    = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
}



VCF_CLASS::VCF_CLASS(){
  etching_version="ETCHING_v1.1.4b (released 2020.12.29)";
}

VCF_CLASS::VCF_CLASS(const std::string infile){
  etching_version="ETCHING_v1.1.4b (released 2020.12.29)";
  read_vcf_file(infile);
}

VCF_CLASS::~VCF_CLASS(){
  clear();
}

void VCF_CLASS::read_vcf_file(const std::string infile){
  std::string tmp;
  VCF vcf;
  std::ifstream fin ( infile );
  while ( std::getline ( fin , tmp ) ){
    if ( tmp[1] == '#' ){
      metainfo += tmp + "\n";
    }
    else{
      header = tmp + "\n";
      break;
    }
  }

  build_id_ref_map(infile);

  
  if ( metainfo.find("##source=ETCHING") != std::string::npos ){
    while ( std::getline ( fin , tmp ) ){
      vcf.parse_etching(tmp);
      insert(vcf);
    }
  }
  else{
    while ( std::getline ( fin , tmp ) ){
      vcf.parse(tmp);
      insert(vcf);
    }
  }
  
  fin.close();
}


void VCF_CLASS::build_id_ref_map(const std::string infile){
  std::ifstream fin ( infile.c_str() );
  std::size_t found;
  std::string tmp;
  std::string id;
  int count = 0;

  const std::string key="##contig=<ID=";

  while ( std::getline ( fin , tmp ) ){
    if ( tmp[1] != '#' ) break;
    found = tmp.find(key);
    if ( found != std::string::npos ){
      tmp = tmp.substr(found+key.size());
      found = tmp.find(",");
      tmp = tmp.substr(0,found);
      id_ref_map[tmp]=count;
      ref_id_map[count]=tmp;
      count ++;
    }
  }
  
  fin.close();
}



void VCF_CLASS::clear(){
  vcf_map.clear();
  genome_info.clear();
  vcf_file.clear();
  reference.clear();
  header.clear();
}


void VCF_CLASS::insert(VCF vcf){
  Position Pos1;
  if ( id_ref_map.find(vcf.chr1) == id_ref_map.end()){
    int Size = id_ref_map.size();
    id_ref_map[vcf.chr1] = Size;
    ref_id_map[Size] = vcf.chr1;
  }
  Pos1.first = id_ref_map[vcf.chr1];
  Pos1.second= vcf.pos1;
  vcf_map[Pos1].push_back(vcf);
}



void VCF_CLASS::insert ( std::string chr1, int pos1, 
			 std::string chr2, int pos2, 
			 std::string sv_id, std::string mate_id, 
			 std::string strand, int sr_val, std::string svtype){
  VCF input;
  Position Pos;

  if ( sr_val < 0 ) sr_val = - sr_val;

  input.chr1 = chr1;
  input.pos1 = pos1;
  input.chr2 = chr2;
  input.pos2 = pos2;
  input.sv_id = sv_id;
  input.mate_id = mate_id;
  input.strand = strand;
  input.sr = sr_val;
  input.svtype = svtype;
  
  input.qual = 30;
  input.filter = ".";

  input.ref = GetNucl(chr1,input.pos1);
  if ( strand == "FR" ){
    input.alt = input.ref + "[" + chr2 + ":" + std::to_string(pos2) + "[";
  }
  else if ( strand == "FF" ){
    input.alt   = input.ref + "]" + chr2 + ":" + std::to_string(pos2) + "]";
  }
  else if ( strand == "RF" ){
    input.alt   = "]" + chr2 + ":" + std::to_string(pos2) + "]" + input.ref ;
  }
  else if ( strand == "RR" ){
    input.alt   = "[" + chr2 + ":" + std::to_string(pos2) + "[" + input.ref ;
  }
  else if ( strand == "F" ) input.alt = input.ref + ".";
  else input.alt = "." + input.ref ;


  insert ( input );
}




VCF_MAP::iterator VCF_CLASS::begin(){return vcf_map.begin();};


VCF_MAP::iterator VCF_CLASS::end(){return vcf_map.end();};

VCF_MAP::iterator VCF_CLASS::find(Position Pos){
  return vcf_map.find(Pos);
}

bool VCF_CLASS::check_vcf( std::string chr1, int pos1,
			   std::string chr2, int pos2,
			   std::string strand ){
  Position Pos1, Pos2;
  Pos1.first  = id_ref_map[chr1];
  Pos1.second = pos1;
  Pos2.first  = id_ref_map[chr2];
  Pos2.second = pos2;
  if ( find(Pos1) != end() ){
    for ( auto & vcf : vcf_map[Pos1] ){
      if ( vcf.chr2 == chr2 && vcf.pos2 == pos2 && vcf.strand == strand ){
	return 1;
      }
    }
  }
  return 0;
}

std::vector < VCF > & VCF_CLASS::operator [](Position Pos){
  return vcf_map[Pos];
}




void VCF_CLASS::get_genome(std::string infile){
  reference = infile;
  get_genome();
}

void VCF_CLASS::get_genome(){
  std::cout << "Reading genome: " << reference << "\n";
  std::ifstream fin(reference.c_str());
  std::string id;
  std::string seq;
  std::string tmp;
  getline ( fin , tmp );
  id = tmp.substr(1);
  while ( getline ( fin , tmp ) ){
    if ( tmp[0] == '>' ){
      genome[id]=seq;
      genome_info.push_back(std::make_pair(id,seq.size()));
      id = tmp.substr(1);
      seq.clear();
    }
    else{
      seq += tmp ;
    }
  }
  genome[id]=seq;
  genome_info.push_back(std::make_pair(id,seq.size()));
  
  fin.close();
}

char VCF_CLASS::GetNucl(const std::string & refChr, const std::size_t & refPos){
  return genome[refChr][refPos-1];
}




void VCF_CLASS::write(){
  std::cout << metainfo << "\n";
  std::cout << header << "\n";
  for ( auto & i : vcf_map ){
    for ( auto & j : i.second ){
      std::cout << j.to_string() << "\n";
    }
  }
}


void VCF_CLASS::fwrite(std::string outfile){
  std::ofstream fout (outfile.c_str());
  fout << metainfo << "\n";
  fout << header << "\n";
  for ( auto & i : vcf_map ){
    for ( auto & j : i.second ){
      fout << j.to_string() << "\n";
    }
  }

  fout.close();
}

void VCF_CLASS::write_short(){
  std::cout << metainfo << "\n";
  std::cout << header << "\n";
  for ( auto & i : vcf_map ){
    for ( auto & j : i.second ){
      std::cout << j.to_string_short() << "\n";
    }
  }
}


void VCF_CLASS::fwrite_short(std::string outfile){
  std::ofstream fout (outfile.c_str());
  fout << metainfo << "\n";
  fout << header << "\n";
  for ( auto & i : vcf_map ){
    for ( auto & j : i.second ){
      fout << j.to_string_short() << "\n";
    }
  }

  fout.close();
}



std::size_t VCF_CLASS::size(){
  std::size_t count = 0 ;
  for ( auto & i : vcf_map ){
    count += i.second.size();
  }
  return count;
}
 
 


void VCF_CLASS::add_features_in_id(){
  for ( auto & i : vcf_map ){
    for ( auto & j : i.second ){
      j.sv_id
        +=":" + std::to_string(j.cr)
	+ ":" + std::to_string(j.sr)
	+ ":" + std::to_string(j.pe)
        + ":" + std::to_string(j.mq)
        + ":" + std::to_string(j.depdif)
        + ":" + std::to_string(j.nxa)
        + ":" + std::to_string(j.tcb)
        + ":" + std::to_string(j.entropy)
        + ":" + std::to_string(j.purity)
        + ":" + std::to_string(j.seqdep);
    }
  }
}






void copy_info ( VCF_CLASS & source, VCF_CLASS & target){
  target = source;
  target.vcf_map.clear();
}




VCF_CLASS typing_SV(VCF_CLASS & input){
  VCF_CLASS output;
  copy_info ( input, output );
  Position Pos1;

  int ins_count = 0;
  int del_count = 0;
  int dup_count = 0;
  int inv_count = 0;

  VCF typed_vcf;
  VCF vcf1;

  std::set < std::string > eliminated;

  for ( auto & j : input.vcf_map ){
    for ( auto & i : j.second ){
      vcf1 = i;

      Pos1.first = input.id_ref_map[vcf1.chr1];
      Pos1.second = vcf1.pos1;

      typed_vcf = vcf1;

      if ( eliminated.find(vcf1.sv_id) == eliminated.end() ){

        if ( vcf1.chr1 == vcf1.chr2 && vcf1.svtype == "BND"){ // intra-Chr. BNDs are classified into DEL, DUP, or INV.
	  if ( vcf1.pos1 == vcf1.pos2 ){
	    if ( vcf1.strand == "FR" ) { // INS
	      ins_count ++ ;
	      
	      typed_vcf.sv_id = "INS" + std::to_string(ins_count);
	      typed_vcf.alt = "<INS>";
	      typed_vcf.svtype = "INS";
	      typed_vcf.mate_id = vcf1.sv_id + "," + vcf1.mate_id;
	      
	      eliminated.insert(vcf1.sv_id);
	      eliminated.insert(vcf1.mate_id);
	    }
	  }
	  else{
	    if ( vcf1.pos1 < vcf1.pos2 ){
	      if ( vcf1.strand == "FR" ) { // DEL

		del_count ++ ;

		typed_vcf.sv_id = "DEL" + std::to_string(del_count);
		typed_vcf.alt = "<DEL>";
		typed_vcf.svtype = "DEL";
		typed_vcf.mate_id = vcf1.sv_id + "," + vcf1.mate_id;

		typed_vcf.sv_id = "DEL" + std::to_string(del_count);

		eliminated.insert(vcf1.sv_id);
		eliminated.insert(vcf1.mate_id);
	      }
	      else if ( vcf1.strand == "RF" ){ // DUP
		dup_count ++ ;

		typed_vcf = vcf1;
		typed_vcf.sv_id = "DUP" + std::to_string(dup_count);
		typed_vcf.alt = "<DUP>";
		typed_vcf.svtype = "DUP";
		typed_vcf.mate_id = vcf1.sv_id + "," + vcf1.mate_id;

		eliminated.insert(vcf1.sv_id);
		eliminated.insert(vcf1.mate_id);
	      }
	      else { // INV
		inv_count ++ ;

		typed_vcf = vcf1;
		typed_vcf.sv_id = "INV" + std::to_string(inv_count);
		typed_vcf.alt = "<INV>";
		typed_vcf.svtype = "INV";
		typed_vcf.mate_id = vcf1.sv_id + "," + vcf1.mate_id;

		eliminated.insert(vcf1.sv_id);
		eliminated.insert(vcf1.mate_id);
	      }
	    }
	    typed_vcf.modify_svtype_info();
	  }
	}
	typed_vcf.make_info();
        output.vcf_map[Pos1].push_back(typed_vcf);
      }
    }
  }
  return output;
}


VCF_CLASS typing_SV_general(VCF_CLASS & input){
  VCF_CLASS output;
  copy_info ( input, output );
  Position Pos1;

  int ins_count = 0;
  int del_count = 0;
  int dup_count = 0;
  int inv_count = 0;

  VCF typed_vcf;
  VCF vcf1;

  std::set < std::string > eliminated;

  for ( auto & j : input.vcf_map ){
    for ( auto & i : j.second ){
      vcf1 = i;

      Pos1.first = input.id_ref_map[vcf1.chr1];
      Pos1.second = vcf1.pos1;

      typed_vcf = vcf1;

      if ( eliminated.find(vcf1.sv_id) == eliminated.end() ){

        if ( vcf1.chr1 == vcf1.chr2 && vcf1.svtype == "BND"){ // intra-Chr. BNDs are classified into DEL, DUP, or INV.
	  if ( vcf1.pos1 == vcf1.pos2 ){
	    if ( vcf1.strand == "FR" ) { // INS // fixed
	      ins_count ++ ;
	      
	      typed_vcf.sv_id = "INS" + std::to_string(ins_count);
	      typed_vcf.alt = "<INS>";
	      typed_vcf.svtype = "INS";
	      typed_vcf.mate_id = vcf1.sv_id + "," + vcf1.mate_id;
	      
	      eliminated.insert(vcf1.sv_id);
	      eliminated.insert(vcf1.mate_id);
	    }
	  }
	  else{
	    if ( vcf1.pos1 < vcf1.pos2 ){
	      if ( vcf1.strand == "FR" ) { // DEL // fixed

		del_count ++ ;

		typed_vcf.sv_id = "DEL" + std::to_string(del_count);
		typed_vcf.alt = "<DEL>";
		typed_vcf.svtype = "DEL";
		typed_vcf.mate_id = vcf1.sv_id + "," + vcf1.mate_id;

		typed_vcf.sv_id = "DEL" + std::to_string(del_count);

		eliminated.insert(vcf1.sv_id);
		eliminated.insert(vcf1.mate_id);
	      }
	      else if ( vcf1.strand == "RF" ){ // DUP //fixed
		dup_count ++ ;

		typed_vcf = vcf1;
		typed_vcf.sv_id = "DUP" + std::to_string(dup_count);
		typed_vcf.alt = "<DUP>";
		typed_vcf.svtype = "DUP";
		typed_vcf.mate_id = vcf1.sv_id + "," + vcf1.mate_id;

		eliminated.insert(vcf1.sv_id);
		eliminated.insert(vcf1.mate_id);
	      }
	      else { // INV
		inv_count ++ ;

		typed_vcf = vcf1;
		typed_vcf.sv_id = "INV" + std::to_string(inv_count);
		typed_vcf.alt = "<INV>";
		typed_vcf.svtype = "INV";
		typed_vcf.mate_id = vcf1.sv_id + "," + vcf1.mate_id;

		eliminated.insert(vcf1.sv_id);
		eliminated.insert(vcf1.mate_id);
	      }
	    }
	    typed_vcf.modify_svtype_info();
	  }
	}
	typed_vcf.make_info();
        output.vcf_map[Pos1].push_back(typed_vcf);
      }
    }
  }
  return output;
}


double return_depdif(Position Pos, std::vector < double > & dep_vec, int Size, double read_length){
  int start = Pos.second - read_length - 1;
  int center = Pos.second - 1;
  int end = Pos.second + read_length - 1;
  
  int max = dep_vec.size();

  if ( start < 0  ) start = 0    ;
  if ( end > Size ) end   = Size ;

  if ( center > max ) center = max;
  if ( end > max ) end = max;

  double left = 0 ;
  double right = 0;
  for ( int pos = start ; pos < center ; pos ++ ) left += dep_vec[pos];
  for ( int pos = center + 1 ; pos < end ; pos ++ ) right += dep_vec[pos];

  double output = ( left - right ) / read_length;
  if ( output < 0 ) output *= -1;
  return output;
}


void VCF_CLASS::calc_features(const std::string input_bam, const int read_length, const int insert_size, const int confi_window){
  
  Position Pos;
  Position Pos1;
  Position Pos2;
  std::string tag;

  BamTools::BamReader reader;
  BamTools::RefVector references;
  BamTools::BamAlignment al;
  BamTools::BamAlignment exal;

  std::map < Position,  int > cr_map;
  std::map < Position, std::map < Position , int > > sr_map;
  std::map < Position, std::map < Position , int > > pe_map;
  std::map < Position,  int > mq_map; //
  std::map < Position,  int > tcb_map; //

  std::map < Position,  int > un_map; // for single break end

  std::vector < double > dep_vec;

  std::map < int , std::set < Position > > ref_pos_map;

  int pe_window_size = 100;

  bool check = false ;
  
  /////////////////////////////////////////////////////////////////
  
  reader.Open(input_bam);
  if ( !reader.OpenIndex(input_bam + ".bai") ) {
    reader.CreateIndex();
  }
  references = reader.GetReferenceData();

  
  for ( auto & i : vcf_map ){
    Pos1 = i.first;
    ref_pos_map[Pos1.first].insert(Pos1);
    for ( auto & j : i.second ){
      j.cr = 0 ;
      // j.sr = 0 ;
      j.pe = 0 ;
      j.mq = 0 ;
      j.nxa = 0;
      j.tcb = 0;
      j.depdif = 0 ;
      j.entropy = 0 ;

      j.cr2 = 0 ;
      // j.sr2 = 0 ;
      j.pe2 = 0 ;
      j.mq2 = 0 ;
      j.nxa2 = 0;
      j.tcb2 = 0;
      j.depdif2 = 0 ;
      j.entropy2 = 0 ;
    }
  }

  /////////////////////////////////////////////////////////////////

  while ( reader.GetNextAlignment(al) ){

    // Update DEPDIF after scanning one chromosome
    if (al.RefID != exal.RefID && check ){
      std::cout << ref_id_map[exal.RefID] <<" done\n";
      int Size = (int) references[exal.RefID].RefLength;
      for ( auto & i : ref_pos_map[exal.RefID] ){
	Pos = i;
	double depdif = return_depdif(Pos, dep_vec, Size, read_length);
	for ( auto & vcf : vcf_map[Pos] ){
	  if ( vcf.svtype != "SND" ){
	    vcf.depdif = depdif ; // update depdif
	    Pos.first = id_ref_map[vcf.chr2];
	    Pos.second= vcf.pos2;
	    vcf.depdif2 = return_depdif(Pos, dep_vec, Size, read_length);
	  }
	  else{
	    vcf.depdif2 = vcf.depdif;
	  }
	}
      }
      dep_vec.clear();
      Size = (int) references[al.RefID].RefLength;
      dep_vec.resize(Size);
    }
    if ( ! check ){
      int Size = (int) references[al.RefID].RefLength;
      dep_vec.resize(Size);
      check = true;
    }
    
    
    /////////////////////////////////////////////////////
    // Scanning clipped read map and mapping quality map
    if ( al.CigarData[0].Type == 'S' || al.CigarData[0].Type == 'H' ){
      Pos1.first = al.RefID;
      Pos1.second = al.Position;
      if ( Pos1.first >= 0 && Pos1.second >= 0 ){
	cr_map [ Pos1 ] ++;
	mq_map [ Pos1 ] += al.MapQuality ;
	tcb_map [ Pos1 ]+= al.CigarData[0].Length ; 
      }
    }
    if ( al.CigarData[al.CigarData.size()-1].Type == 'S' || al.CigarData[al.CigarData.size()-1].Type == 'H' ){
      Pos1.first = al.RefID;
      Pos1.second = al.GetEndPosition() - 1;
      if ( Pos1.first >= 0 && Pos1.second >= 0 ){
	cr_map [ Pos1 ] ++ ;
	mq_map [ Pos1 ] += al.MapQuality ;
	tcb_map [ Pos1 ]+= al.CigarData[al.CigarData.size()-1].Length ; 
      }
    }

    /*
    /////////////////////////////////////////////////////
    // Scanning split read map
    std::string Tag;
    if ( al.GetTag<std::string>("SA",Tag) ){
      std::vector < std::string > tag_vec;
      while ( Tag.size() > 0 ){
	std::string tmp = parse_SA_tag(Tag);
	tag_vec.push_back( tmp );
      }
      for ( auto tag : tag_vec){
	if ( al.CigarData[al.CigarData.size()-1].Type == 'S' || al.CigarData[al.CigarData.size()-1].Type == 'H' ){
	  Pos1.first = al.RefID ;
	  Pos1.second = al.GetEndPosition()-1;
	  Pos2 = find_split(tag,id_ref_map);
	  if ( Pos2.second > -1 ){
	    sr_map[Pos1][Pos2]++;
	  }
	}
	if ( al.CigarData[0].Type == 'S' || al.CigarData[0].Type == 'H' ){
	  Pos1.first = al.RefID ;
	  Pos1.second = al.Position;
	  Pos2 = find_split(tag,id_ref_map);
	  if ( Pos2.second > -1 ){
	    sr_map[Pos1][Pos2]++;
	  }
	}
      }
    }
    
    // SND
    if ( ! al.IsMateMapped () ){
      if ( al.CigarData[al.CigarData.size()-1].Type == 'S' || al.CigarData[al.CigarData.size()-1].Type == 'H' ){
	Pos1.first = al.RefID ;
	Pos1.second = al.GetEndPosition()-1;
	un_map[Pos1] ++;
      }
      if ( al.CigarData[0].Type == 'S' || al.CigarData[0].Type == 'H' ){
	Pos1.first = al.RefID ;
	Pos1.second = al.Position;
	un_map[Pos1] ++;
      }
    }
    */

    /////////////////////////////////////////////////////
    // Scanning discordant map
    Pos1.first = al.RefID;
    Pos1.second = al.Position + read_length / 2 ;
    Pos1.second /= pe_window_size ;

    Pos2.first = al.MateRefID;
    Pos2.second = al.MatePosition + read_length / 2 ;
    Pos2.second /= pe_window_size ;

    if ( Pos1.first != Pos2.first ){
      pe_map[Pos1][Pos2]=1;
      pe_map[Pos2][Pos1]=1;
      // pe_map[Pos1][Pos2]++;
      // pe_map[Pos2][Pos1]++;
    }
    else if ( Pos1.second - Pos2.second > 2 * insert_size / pe_window_size || Pos1.second - Pos2.second < - 2 * insert_size / pe_window_size ){
      pe_map[Pos1][Pos2]=1;
      pe_map[Pos2][Pos1]=1;
      // pe_map[Pos1][Pos2]++;
      // pe_map[Pos2][Pos1]++;
    }



    /////////////////////////////////////////////////////
    // depth array
    const int start = al.Position ;
    const int end = al.GetEndPosition() - 1 ;

    for ( int i = start ; i < end ; i ++ ){
      dep_vec[i] ++ ;
    }

    exal = al;
  } // End of while()

  // Update DEPDIF after scanning one chromosome
  int Size = (int) references[exal.RefID].RefLength;
  for ( auto & i : ref_pos_map[exal.RefID] ){
    Pos = i;
    double depdif = return_depdif(Pos, dep_vec, Size, read_length);
    for ( auto & vcf : vcf_map[Pos] ){
      vcf.depdif = depdif ; // update depdif
      Pos.first = id_ref_map[vcf.chr2];
      Pos.second= vcf.pos2;
      vcf.depdif2 = return_depdif(Pos, dep_vec, Size, read_length);
    }
  }
  dep_vec.clear();

  std::cout << "Scanning done\n";
  std::cout << "Complete: DEPDIF\n";

  
  /////////////////////////////////////////////////////
  // Calculatiing CR
  for ( auto & i : vcf_map ){
    // cr
    Pos1 = i.first;
    Position Posa = Pos1;
    for ( auto & vcf : i.second ){
      for (int k = - confi_window ; k < confi_window - 1; k ++ ){
	Posa.second = Pos1.second + k ;
	if ( cr_map.find(Posa) != cr_map.end() ){
	  vcf.cr += cr_map [Posa];
	}
	if ( mq_map.find(Posa) != mq_map.end() ){
	  vcf.mq += mq_map [Posa];
	}
	if ( tcb_map.find(Posa) != tcb_map.end() ){
	  vcf.tcb+= tcb_map[Posa];
	}
      }
      if ( vcf.svtype != "SND" ){
	Pos2.first = id_ref_map[vcf.chr2];
	Pos2.second = vcf.pos2;
	Position Posb = Pos2;
	for (int k = - confi_window ; k < confi_window - 1; k ++ ){
	  Posb.second = Pos2.second + k ;
	  if ( cr_map.find(Posb) != cr_map.end() ){
	    vcf.cr2 += cr_map [Posb];
	  }
	  if ( mq_map.find(Posb) != mq_map.end() ){
	    vcf.mq2 += mq_map [Posb];
	  }
	  if ( tcb_map.find(Posb) != tcb_map.end() ){
	    vcf.tcb2+= tcb_map[Posb];
	  }
	}
      }
      else{
	vcf.cr2 = vcf.cr;
	vcf.mq2 = vcf.mq;
	vcf.tcb2= vcf.tcb;
      }
    }
  }


  /////////////////////////////////////////////////////
  // Calculatiing MQ
  for ( auto & i : vcf_map ){
    for ( auto & vcf : i.second ){
      if ( vcf.cr !=0 ) vcf.mq  /= vcf.cr ;
      if ( vcf.cr2!=0 ) vcf.mq2 /= vcf.cr2;
    }
  }
  std::cout << "Complete: CR\n";
  std::cout << "Complete: MQ\n";
  std::cout << "Complete: TCB\n";

  /*
  /////////////////////////////////////////////////////
  // Calculatiing SR
  for ( auto & i : vcf_map ){
    Pos1 = i.first;
    Position Posa = Pos1;
    for ( auto & vcf : i.second ){
      if ( vcf.svtype!="SND"){
	Pos2.first = id_ref_map[vcf.chr2];
	Pos2.second = vcf.pos2;
	Position Posb = Pos2;
	for (int k = - confi_window ; k < confi_window - 1; k ++ ){
	  Posa.second = Pos1.second + k ;
	  if ( sr_map.find(Posa) != sr_map.end() ){
	    for ( int l = - confi_window ; l <= confi_window - 1; l ++ ){
	      Posb.second = Pos2.second + l ;
	      if ( sr_map[Posa].find(Posb) != sr_map[Posa].end() ){
		vcf.sr += sr_map[Posa][Posb];
	      }
	    }
	  }
	}
      }
      else{ // SND
	for (int k = - confi_window ; k < confi_window - 1; k ++ ){ 
	  Posa.second = Pos1.second + k ;
	  if ( un_map.find(Posa) != un_map.end() ){ 
	    vcf.sr += un_map[Posa];
	    vcf.pe += un_map[Posa];
	  }
	}
      }
    }
  }
  */
  for ( auto & i : vcf_map ){
    for ( auto & vcf : i.second ){
      vcf.sr2 = vcf.sr;
    }
  }

  std::cout << "Complete: SR\n";


  /////////////////////////////////////////////////////
  // Calculatiing PE
  for ( auto & i : vcf_map ){
    Pos1 = i.first;
    Position Posa = Pos1;
    Posa.second += read_length / 2;
    Posa.second /= pe_window_size;
    for ( auto & vcf : i.second ){
      Pos2.first = id_ref_map[vcf.chr2];
      Pos2.second = vcf.pos2;
      Position Posb = Pos2;
      Posb.second += read_length / 2;
      Posb.second /= pe_window_size ;
      for (int k = - 2 * insert_size / pe_window_size ; k <= 2 * insert_size / pe_window_size ; k ++ ){
	Posa.second = ( Pos1.second / pe_window_size ) + k ;
	if ( pe_map.find(Posa) != pe_map.end() ){
	  for ( int l = - 2 * insert_size / pe_window_size ; l <= 2 * insert_size / pe_window_size ; l ++ ){
	    Posb.second = ( Pos2.second / pe_window_size ) + l ;
	    if ( pe_map[Posa].find(Posb) != pe_map[Posa].end() ){
	      vcf.pe += pe_map[Posa][Posb];
	    }
	  }
	}
      }
    }
  }
  for ( auto & i : vcf_map ){
    for ( auto & vcf : i.second ){
      vcf.pe /= 2;
    }
  }

  // If sv type is SND, it is assumed that its PE is same with SR.
  for ( auto & i : vcf_map ){
    for ( auto & vcf : i.second ){
      if ( vcf.svtype == "SND" ){
	vcf.pe = vcf.sr;
      }
    }
  }

  for ( auto & i : vcf_map ){
    for ( auto & vcf : i.second ){
      vcf.pe2 = vcf.pe;
    }
  }

  std::cout << "Complete: PE \n";


  /////////////////////////////////////////////////////
  // Calculatiing ENTROPY
  std::map<Position,double> entropy_map;

  for ( auto & i : vcf_map ){
    Pos1 = i.first;
    for ( auto & vcf : i.second ){
      Pos2.first = id_ref_map[vcf.chr2];
      Pos2.second = vcf.pos2;
      sr_map[Pos1][Pos2]=vcf.sr;
    }
  }

  for ( auto & i : sr_map ){
    if ( i.second.size() < 2 ){
      entropy_map[i.first] = 0 ;
    }
    else{
      double sum = 0, prob = 0, entropy = 0 ;
      for ( auto & j : i.second ){
        sum += (double) j.second ;
      }
      for ( auto & j : i.second ){
        if ( sum > 0 ) prob = ( (double) j.second ) / sum ;
	if ( prob > 0 ) entropy = - prob * log ( prob );
      }
      if ( entropy < 0.000001 ) entropy = 0;
      entropy_map[i.first]=entropy;
    }
  }


  for ( auto & i : vcf_map ){
    Pos1 = i.first;
    Position Posa = Pos1;
    for ( auto & vcf : i.second ){
      for (int k = - confi_window ; k < confi_window - 1; k ++ ){
	Posa.second = Pos1.second + k ;
	if ( entropy_map.find(Posa) != entropy_map.end() ){
	  vcf.entropy += entropy_map[Posa];
	}
      }
      if ( vcf.svtype != "SND" ){
	Pos2.first = id_ref_map[vcf.chr2];
	Pos2.second = vcf.pos2;
	Position Posb = Pos2;
	for (int k = - confi_window ; k < confi_window - 1; k ++ ){
	  Posb.second = Pos2.second + k ;
	  if ( entropy_map.find(Posb) != entropy_map.end() ){
	    vcf.entropy2 += entropy_map[Posb];
	  }
	}
      }
      else{
	vcf.entropy2 = vcf.entropy;
      }
    }
  }
  std::cout << "Complete: ENTROPY\n";



  /////////////////////////////////////////////////////
  // make info field
  for ( auto & i : vcf_map ){
    for ( auto & j : i.second ){
      j.make_info(); // update info
    }
  }
  std::cout << "Complete: feature calculation\n";

}






void VCF_CLASS::calc_features_general(const std::string input_bam, const int read_length, const int insert_size, const int confi_window){
  
  Position Pos;
  Position Pos1;
  Position Pos2;
  std::string tag;

  BamTools::BamReader reader;
  BamTools::RefVector references;
  BamTools::BamAlignment al;
  BamTools::BamAlignment exal;

  std::map < Position,  int > cr_map;
  std::map < Position, std::map < Position , int > > sr_map;
  std::map < Position, std::map < Position , int > > pe_map;
  std::map < Position,  int > mq_map; //
  std::map < Position,  int > tcb_map; //

  std::map < Position,  int > un_map; // for single break end

  std::vector < double > dep_vec;

  std::map < int , std::set < Position > > ref_pos_map;

  int pe_window_size = 100;

  bool check = false ;
  
  /////////////////////////////////////////////////////////////////
  
  reader.Open(input_bam);
  if ( !reader.OpenIndex(input_bam + ".bai") ) {
    reader.CreateIndex();
  }
  references = reader.GetReferenceData();

  
  for ( auto & i : vcf_map ){
    Pos1 = i.first;
    ref_pos_map[Pos1.first].insert(Pos1);
    for ( auto & j : i.second ){
      j.cr = 0 ;
      j.sr = 0 ;
      j.pe = 0 ;
      j.mq = 0 ;
      j.nxa = 0;
      j.tcb = 0;
      j.depdif = 0 ;
      j.entropy = 0 ;

      j.cr2 = 0 ;
      j.sr2 = 0 ;
      j.pe2 = 0 ;
      j.mq2 = 0 ;
      j.nxa2 = 0;
      j.tcb2 = 0;
      j.depdif2 = 0 ;
      j.entropy2 = 0 ;
    }
  }

  /////////////////////////////////////////////////////////////////

  while ( reader.GetNextAlignment(al) ){

    // Update DEPDIF after scanning one chromosome
    if (al.RefID != exal.RefID && check ){
      std::cout << ref_id_map[exal.RefID] <<" done\n";
      int Size = (int) references[exal.RefID].RefLength;
      for ( auto & i : ref_pos_map[exal.RefID] ){
	Pos = i;
	double depdif = return_depdif(Pos, dep_vec, Size, read_length);
	for ( auto & vcf : vcf_map[Pos] ){
	  if ( vcf.svtype != "SND" ){
	    vcf.depdif = depdif ; // update depdif
	    Pos.first = id_ref_map[vcf.chr2];
	    Pos.second= vcf.pos2;
	    vcf.depdif2 = return_depdif(Pos, dep_vec, Size, read_length);
	  }
	  else{
	    vcf.depdif2 = vcf.depdif;
	  }
	}
      }
      dep_vec.clear();
      Size = (int) references[al.RefID].RefLength;
      dep_vec.resize(Size);
    }
    if ( ! check ){
      int Size = (int) references[al.RefID].RefLength;
      dep_vec.resize(Size);
      check = true;
    }
    
    
    /////////////////////////////////////////////////////
    // Scanning clipped read map and mapping quality map
    if ( al.CigarData[0].Type == 'S' || al.CigarData[0].Type == 'H' ){
      Pos1.first = al.RefID;
      Pos1.second = al.Position;
      if ( Pos1.first >= 0 && Pos1.second >= 0 ){
	cr_map [ Pos1 ] ++;
	mq_map [ Pos1 ] += al.MapQuality ;
	tcb_map [ Pos1 ]+= al.CigarData[0].Length ; 
      }
    }
    if ( al.CigarData[al.CigarData.size()-1].Type == 'S' || al.CigarData[al.CigarData.size()-1].Type == 'H' ){
      Pos1.first = al.RefID;
      Pos1.second = al.GetEndPosition() - 1;
      if ( Pos1.first >= 0 && Pos1.second >= 0 ){
	cr_map [ Pos1 ] ++ ;
	mq_map [ Pos1 ] += al.MapQuality ;
	tcb_map [ Pos1 ]+= al.CigarData[al.CigarData.size()-1].Length ; 
      }
    }


    /////////////////////////////////////////////////////
    // Scanning split read map
    std::string Tag;
    if ( al.GetTag<std::string>("SA",Tag) ){
      std::vector < std::string > tag_vec;
      while ( Tag.size() > 0 ){
	std::string tmp = parse_SA_tag(Tag);
	tag_vec.push_back( tmp );
      }
      for ( auto tag : tag_vec){
	if ( al.CigarData[al.CigarData.size()-1].Type == 'S' || al.CigarData[al.CigarData.size()-1].Type == 'H' ){
	  Pos1.first = al.RefID ;
	  Pos1.second = al.GetEndPosition()-1;
	  Pos2 = find_split(tag,id_ref_map);
	  if ( Pos2.second > -1 ){
	    sr_map[Pos1][Pos2]++;
	  }
	}
	if ( al.CigarData[0].Type == 'S' || al.CigarData[0].Type == 'H' ){
	  Pos1.first = al.RefID ;
	  Pos1.second = al.Position;
	  Pos2 = find_split(tag,id_ref_map);
	  if ( Pos2.second > -1 ){
	    sr_map[Pos1][Pos2]++;
	  }
	}
      }
    }
    // SND
    if ( ! al.IsMateMapped () ){
      if ( al.CigarData[al.CigarData.size()-1].Type == 'S' || al.CigarData[al.CigarData.size()-1].Type == 'H' ){
	Pos1.first = al.RefID ;
	Pos1.second = al.GetEndPosition()-1;
	un_map[Pos1] ++;
      }
      if ( al.CigarData[0].Type == 'S' || al.CigarData[0].Type == 'H' ){
	Pos1.first = al.RefID ;
	Pos1.second = al.Position;
	un_map[Pos1] ++;
      }
    }


    /////////////////////////////////////////////////////
    // Scanning discordant map
    Pos1.first = al.RefID;
    Pos1.second = al.Position + read_length / 2 ;
    Pos1.second /= pe_window_size ;

    Pos2.first = al.MateRefID;
    Pos2.second = al.MatePosition + read_length / 2 ;
    Pos2.second /= pe_window_size ;

    if ( Pos1.first != Pos2.first ){
      pe_map[Pos1][Pos2]=1;
      pe_map[Pos2][Pos1]=1;
      // pe_map[Pos1][Pos2]++;
      // pe_map[Pos2][Pos1]++;
    }
    else if ( Pos1.second - Pos2.second > 2 * insert_size / pe_window_size || Pos1.second - Pos2.second < - 2 * insert_size / pe_window_size ){
      pe_map[Pos1][Pos2]=1;
      pe_map[Pos2][Pos1]=1;
      // pe_map[Pos1][Pos2]++;
      // pe_map[Pos2][Pos1]++;
    }



    /////////////////////////////////////////////////////
    // depth array
    const int start = al.Position ;
    const int end = al.GetEndPosition() - 1 ;

    for ( int i = start ; i < end ; i ++ ){
      dep_vec[i] ++ ;
    }

    exal = al;
  } // End of while()

  // Update DEPDIF after scanning one chromosome
  int Size = (int) references[exal.RefID].RefLength;
  for ( auto & i : ref_pos_map[exal.RefID] ){
    Pos = i;
    double depdif = return_depdif(Pos, dep_vec, Size, read_length);
    for ( auto & vcf : vcf_map[Pos] ){
      vcf.depdif = depdif ; // update depdif
      Pos.first = id_ref_map[vcf.chr2];
      Pos.second= vcf.pos2;
      vcf.depdif2 = return_depdif(Pos, dep_vec, Size, read_length);
    }
  }
  dep_vec.clear();

  std::cout << "Scanning done\n";
  std::cout << "Complete: DEPDIF\n";

  
  /////////////////////////////////////////////////////
  // Calculatiing CR
  for ( auto & i : vcf_map ){
    // cr
    Pos1 = i.first;
    Position Posa = Pos1;
    for ( auto & vcf : i.second ){
      for (int k = - confi_window ; k < confi_window - 1; k ++ ){
	Posa.second = Pos1.second + k ;
	if ( cr_map.find(Posa) != cr_map.end() ){
	  vcf.cr += cr_map [Posa];
	  vcf.mq += mq_map [Posa];
	  vcf.tcb+= tcb_map[Posa];
	}
      }
      if ( vcf.svtype != "SND" ){
	Pos2.first = id_ref_map[vcf.chr2];
	Pos2.second = vcf.pos2;
	Position Posb = Pos2;
	for (int k = - confi_window ; k < confi_window - 1; k ++ ){
	  Posb.second = Pos2.second + k ;
	  if ( cr_map.find(Posb) != cr_map.end() ){
	    vcf.cr2 += cr_map [Posb];
	    vcf.mq2 += mq_map [Posb];
	    vcf.tcb2+= tcb_map[Posb];
	  }
	}
      }
      else{
	vcf.cr2 = vcf.cr;
	vcf.mq2 = vcf.mq;
	vcf.tcb2= vcf.tcb;
      }
    }
  }


  /////////////////////////////////////////////////////
  // Calculatiing MQ
  for ( auto & i : vcf_map ){
    for ( auto & vcf : i.second ){
      if ( vcf.cr !=0 ) vcf.mq  /= vcf.cr ;
      if ( vcf.cr2!=0 ) vcf.mq2 /= vcf.cr2;
    }
  }
  std::cout << "Complete: CR\n";
  std::cout << "Complete: MQ\n";
  std::cout << "Complete: TCB\n";


  /////////////////////////////////////////////////////
  // Calculatiing SR
  for ( auto & i : vcf_map ){
    Pos1 = i.first;
    Position Posa = Pos1;
    for ( auto & vcf : i.second ){
      if ( vcf.svtype!="SND"){
	Pos2.first = id_ref_map[vcf.chr2];
	Pos2.second = vcf.pos2;
	Position Posb = Pos2;
	for (int k = - confi_window ; k < confi_window - 1; k ++ ){
	  Posa.second = Pos1.second + k ;
	  if ( sr_map.find(Posa) != sr_map.end() ){
	    for ( int l = - confi_window ; l <= confi_window - 1; l ++ ){
	      Posb.second = Pos2.second + l ;
	      if ( sr_map[Posa].find(Posb) != sr_map[Posa].end() ){
		vcf.sr += sr_map[Posa][Posb];
	      }
	    }
	  }
	}
      }
      else{ // SND
	for (int k = - confi_window ; k < confi_window - 1; k ++ ){ 
	  Posa.second = Pos1.second + k ;
	  if ( un_map.find(Posa) != un_map.end() ){ 
	    vcf.sr += un_map[Posa];
	    vcf.pe += un_map[Posa];
	  }
	}
      }
    }
  }

  for ( auto & i : vcf_map ){
    for ( auto & vcf : i.second ){
      vcf.sr2 = vcf.sr;
    }
  }

  std::cout << "Complete: SR\n";


  /////////////////////////////////////////////////////
  // Calculatiing PE
  for ( auto & i : vcf_map ){
    Pos1 = i.first;
    Position Posa = Pos1;
    Posa.second += read_length / 2;
    Posa.second /= pe_window_size;
    for ( auto & vcf : i.second ){
      Pos2.first = id_ref_map[vcf.chr2];
      Pos2.second = vcf.pos2;
      Position Posb = Pos2;
      Posb.second += read_length / 2;
      Posb.second /= pe_window_size ;
      for (int k = - 2 * insert_size / pe_window_size ; k <= 2 * insert_size / pe_window_size ; k ++ ){
	Posa.second = ( Pos1.second / pe_window_size ) + k ;
	if ( pe_map.find(Posa) != pe_map.end() ){
	  for ( int l = - 2 * insert_size / pe_window_size ; l <= 2 * insert_size / pe_window_size ; l ++ ){
	    Posb.second = ( Pos2.second / pe_window_size ) + l ;
	    if ( pe_map[Posa].find(Posb) != pe_map[Posa].end() ){
	      vcf.pe += pe_map[Posa][Posb];
	    }
	  }
	}
      }
    }
  }
  for ( auto & i : vcf_map ){
    for ( auto & vcf : i.second ){
      vcf.pe /= 2;
    }
  }

  // If sv type is SND, it is assumed that its PE is same with SR.
  for ( auto & i : vcf_map ){
    for ( auto & vcf : i.second ){
      if ( vcf.svtype == "SND" ){
	vcf.pe = vcf.sr;
      }
    }
  }

  for ( auto & i : vcf_map ){
    for ( auto & vcf : i.second ){
      vcf.pe2 = vcf.pe;
    }
  }

  std::cout << "Complete: PE \n";


  /////////////////////////////////////////////////////
  // Calculatiing ENTROPY
  std::map<Position,double> entropy_map;
  for ( auto & i : sr_map ){
    if ( i.second.size() < 2 ){
      entropy_map[i.first] = 0 ;
    }
    else{
      double sum = 0, prob = 0, entropy = 0 ;
      for ( auto & j : i.second ){
        sum += (double) j.second ;
      }
      for ( auto & j : i.second ){
        if ( sum > 0 ) prob = ( (double) j.second ) / sum ;
	if ( prob > 0 ) entropy = - prob * log ( prob );
      }
      if ( entropy < 0.000001 ) entropy = 0;
      entropy_map[i.first]=entropy;
    }
  }


  for ( auto & i : vcf_map ){
    Pos1 = i.first;
    Position Posa = Pos1;
    for ( auto & vcf : i.second ){
      for (int k = - confi_window ; k < confi_window - 1; k ++ ){
	Posa.second = Pos1.second + k ;
	if ( entropy_map.find(Posa) != entropy_map.end() ){
	  vcf.entropy += entropy_map[Posa];
	}
      }
      if ( vcf.svtype != "SND" ){
	Pos2.first = id_ref_map[vcf.chr2];
	Pos2.second = vcf.pos2;
	Position Posb = Pos2;
	for (int k = - confi_window ; k < confi_window - 1; k ++ ){
	  Posb.second = Pos2.second + k ;
	  if ( entropy_map.find(Posb) != entropy_map.end() ){
	    vcf.entropy2 += entropy_map[Posb];
	  }
	}
      }
      else{
	vcf.entropy2 = vcf.entropy;
      }
    }
  }
  std::cout << "Complete: ENTROPY\n";



  /////////////////////////////////////////////////////
  // make info field
  for ( auto & i : vcf_map ){
    for ( auto & j : i.second ){
      j.make_info(); // update info
    }
  }
  std::cout << "Complete: feature calculation\n";

}






void VCF_CLASS::make_info(){
  for ( auto & i : vcf_map ){
    for ( auto & j : i.second ){
      j.make_info();
    }
  }  
}

void calc_features(const std::string input_bam, std::vector < VCF_CLASS > & container_vec,
		   const int read_length, const int insert_size, const int confi_window){
  
  Position Pos;
  Position Pos1;
  Position Pos2;
  std::string tag;

  BamTools::BamReader reader;
  BamTools::RefVector references;
  BamTools::BamAlignment al;
  BamTools::BamAlignment exal;

  std::map < Position,  int > cr_map;
  std::map < Position, std::map < Position , int > > sr_map;
  std::map < Position, std::map < Position , int > > pe_map;
  std::map < Position,  int > mq_map; //
  std::map < Position,  int > tcb_map; //

  std::map < Position,  int > un_map; // for single break end

  std::vector < double > dep_vec;

  // std::map < int , std::set < Position > > ref_pos_map;
  std::vector < std::map < int , std::set < Position > > > ref_pos_map ( container_vec.size() );
  
  int pe_window_size = 100;

  bool check = false ;
  

  std::unordered_map < std::string , int > id_ref_map = container_vec[0].id_ref_map;
  std::unordered_map < int , std::string > ref_id_map = container_vec[0].ref_id_map;
  /////////////////////////////////////////////////////////////////
  
  reader.Open(input_bam);
  if ( !reader.OpenIndex(input_bam + ".bai") ) {
    reader.CreateIndex();
  }
  references = reader.GetReferenceData();

  for ( std::size_t con = 0 ; con < container_vec.size() ; con ++ ){
    for ( auto & i : container_vec[con].vcf_map ){
      Pos1 = i.first;
      ref_pos_map[con][Pos1.first].insert(Pos1);
      for ( auto & j : i.second ){
	j.cr = 0 ;
	j.sr = 0 ;
	j.pe = 0 ;
	j.mq = 0 ;
	j.nxa = 0;
	j.tcb = 0;
	j.depdif = 0 ;
	j.entropy = 0 ;

	j.cr2 = 0 ;
	j.sr2 = 0 ;
	j.pe2 = 0 ;
	j.mq2 = 0 ;
	j.nxa2 = 0;
	j.tcb2 = 0;
	j.depdif2 = 0 ;
	j.entropy2 = 0 ;
      }
    }
  }

  /////////////////////////////////////////////////////////////////

  while ( reader.GetNextAlignment(al) ){

    // Update DEPDIF after scanning one chromosome
    if (al.RefID != exal.RefID && check ){
      std::cout << ref_id_map[exal.RefID] <<" done\n";
      int Size = (int) references[exal.RefID].RefLength;

      for ( std::size_t con = 0 ; con < container_vec.size() ; con ++ ){
	for ( auto & i : ref_pos_map[con][exal.RefID] ){
	  Pos = i;
	  double depdif = return_depdif(Pos, dep_vec, Size, read_length);
	  for ( auto & vcf : container_vec[con].vcf_map[Pos] ){
	    if ( vcf.svtype != "SND" ){
	      vcf.depdif = depdif ; // update depdif
	      Pos.first = id_ref_map[vcf.chr2];
	      Pos.second= vcf.pos2;
	      vcf.depdif2 = return_depdif(Pos, dep_vec, Size, read_length);
	    }
	    else{
	      vcf.depdif2 = vcf.depdif;
	    }
	  }
	}
      }
      dep_vec.clear();
      Size = (int) references[al.RefID].RefLength;
      dep_vec.resize(Size);
    }
    if ( ! check ){
      int Size = (int) references[al.RefID].RefLength;
      dep_vec.resize(Size);
      check = true;
    }
    
    
    /////////////////////////////////////////////////////
    // Scanning clipped read map and mapping quality map
    if ( al.CigarData[0].Type == 'S' || al.CigarData[0].Type == 'H' ){
      Pos1.first = al.RefID;
      Pos1.second = al.Position;
      if ( Pos1.first >= 0 && Pos1.second >= 0 ){
	cr_map [ Pos1 ] ++;
	mq_map [ Pos1 ] += al.MapQuality ;
	tcb_map [ Pos1 ]+= al.CigarData[0].Length ; 
      }
    }
    if ( al.CigarData[al.CigarData.size()-1].Type == 'S' || al.CigarData[al.CigarData.size()-1].Type == 'H' ){
      Pos1.first = al.RefID;
      Pos1.second = al.GetEndPosition() - 1;
      if ( Pos1.first >= 0 && Pos1.second >= 0 ){
	cr_map [ Pos1 ] ++ ;
	mq_map [ Pos1 ] += al.MapQuality ;
	tcb_map [ Pos1 ]+= al.CigarData[al.CigarData.size()-1].Length ; 
      }
    }


    /////////////////////////////////////////////////////
    // Scanning split read map
    std::string Tag;
    if ( al.GetTag<std::string>("SA",Tag) ){
      std::vector < std::string > tag_vec;
      while ( Tag.size() > 0 ){
	std::string tmp = parse_SA_tag(Tag);
	tag_vec.push_back( tmp );
      }
      for ( auto tag : tag_vec){
	if ( al.CigarData[al.CigarData.size()-1].Type == 'S' || al.CigarData[al.CigarData.size()-1].Type == 'H' ){
	  Pos1.first = al.RefID ;
	  Pos1.second = al.GetEndPosition()-1;
	  Pos2 = find_split(tag,id_ref_map);
	  if ( Pos2.second > -1 ){
	    sr_map[Pos1][Pos2]++;
	  }
	}
	if ( al.CigarData[0].Type == 'S' || al.CigarData[0].Type == 'H' ){
	  Pos1.first = al.RefID ;
	  Pos1.second = al.Position;
	  Pos2 = find_split(tag,id_ref_map);
	  if ( Pos2.second > -1 ){
	    sr_map[Pos1][Pos2]++;
	  }
	}
      }
    }
    
    // SND
    if ( ! al.IsMateMapped () ){
      if ( al.CigarData[al.CigarData.size()-1].Type == 'S' || al.CigarData[al.CigarData.size()-1].Type == 'H' ){
	Pos1.first = al.RefID ;
	Pos1.second = al.GetEndPosition()-1;
	un_map[Pos1] ++;
      }
      if ( al.CigarData[0].Type == 'S' || al.CigarData[0].Type == 'H' ){
	Pos1.first = al.RefID ;
	Pos1.second = al.Position;
	un_map[Pos1] ++;
      }
    }


    /////////////////////////////////////////////////////
    // Scanning discordant map
    Pos1.first = al.RefID;
    Pos1.second = al.Position + read_length / 2 ;
    Pos1.second /= pe_window_size ;

    Pos2.first = al.MateRefID;
    Pos2.second = al.MatePosition + read_length / 2 ;
    Pos2.second /= pe_window_size ;

    if ( Pos1.first != Pos2.first ){
      pe_map[Pos1][Pos2]=1;
      pe_map[Pos2][Pos1]=1;
      // pe_map[Pos1][Pos2]++;
      // pe_map[Pos2][Pos1]++;
    }
    else if ( Pos1.second - Pos2.second > 2 * insert_size / pe_window_size || Pos1.second - Pos2.second < - 2 * insert_size / pe_window_size ){
      pe_map[Pos1][Pos2]=1;
      pe_map[Pos2][Pos1]=1;
      // pe_map[Pos1][Pos2]++;
      // pe_map[Pos2][Pos1]++;
    }



    /////////////////////////////////////////////////////
    // depth array
    const int start = al.Position ;
    const int end = al.GetEndPosition() - 1 ;

    for ( int i = start ; i < end ; i ++ ){
      dep_vec[i] ++ ;
    }

    exal = al;
  } // End of while()


  // Update DEPDIF after scanning one chromosome
  int Size = (int) references[exal.RefID].RefLength;
  for ( std::size_t con = 0 ; con < container_vec.size() ; con ++ ){
    for ( auto & i : ref_pos_map[con][exal.RefID] ){
      Pos = i;
      double depdif = return_depdif(Pos, dep_vec, Size, read_length);
      for ( auto & vcf : container_vec[con].vcf_map[Pos] ){
	vcf.depdif = depdif ; // update depdif
	Pos.first = id_ref_map[vcf.chr2];
	Pos.second= vcf.pos2;
	vcf.depdif2 = return_depdif(Pos, dep_vec, Size, read_length);
      }
    }
  }
  dep_vec.clear();

  std::cout << "Scanning done\n";
  std::cout << "Complete: DEPDIF\n";

  
  /////////////////////////////////////////////////////
  // Calculatiing CR
  for ( std::size_t con = 0 ; con < container_vec.size() ; con ++ ){
    for ( auto & i : container_vec[con].vcf_map ){
      // cr
      Pos1 = i.first;
      Position Posa = Pos1;
      for ( auto & vcf : i.second ){
	for (int k = - confi_window ; k < confi_window - 1; k ++ ){
	  Posa.second = Pos1.second + k ;
	  if ( cr_map.find(Posa) != cr_map.end() ){
	    vcf.cr += cr_map [Posa];
	    vcf.mq += mq_map [Posa];
	    vcf.tcb+= tcb_map[Posa];
	  }
	}
	if ( vcf.svtype != "SND" ){
	  Pos2.first = id_ref_map[vcf.chr2];
	  Pos2.second = vcf.pos2;
	  Position Posb = Pos2;
	  for (int k = - confi_window ; k < confi_window - 1; k ++ ){
	    Posb.second = Pos2.second + k ;
	    if ( cr_map.find(Posb) != cr_map.end() ){
	      vcf.cr2 += cr_map [Posb];
	      vcf.mq2 += mq_map [Posb];
	      vcf.tcb2+= tcb_map[Posb];
	    }
	  }
	}
	else{
	  vcf.cr2 = vcf.cr;
	  vcf.mq2 = vcf.mq;
	  vcf.tcb2= vcf.tcb;
	}
      }
    }
  }

  /////////////////////////////////////////////////////
  // Calculatiing MQ
  for ( std::size_t con = 0 ; con < container_vec.size() ; con ++ ){
    for ( auto & i : container_vec[con].vcf_map ){
      for ( auto & vcf : i.second ){
	if ( vcf.cr !=0 ) vcf.mq  /= vcf.cr ;
	if ( vcf.cr2!=0 ) vcf.mq2 /= vcf.cr2;
      }
    }
  }
  std::cout << "Complete: CR\n";
  std::cout << "Complete: MQ\n";
  std::cout << "Complete: TCB\n";


  /////////////////////////////////////////////////////
  // Calculatiing SR
  for ( std::size_t con = 0 ; con < container_vec.size() ; con ++ ){
    for ( auto & i : container_vec[con].vcf_map ){
      Pos1 = i.first;
      Position Posa = Pos1;
      for ( auto & vcf : i.second ){
	if ( vcf.svtype!="SND"){
	  Pos2.first = id_ref_map[vcf.chr2];
	  Pos2.second = vcf.pos2;
	  Position Posb = Pos2;
	  for (int k = - confi_window ; k < confi_window - 1; k ++ ){
	    Posa.second = Pos1.second + k ;
	    if ( sr_map.find(Posa) != sr_map.end() ){
	      for ( int l = - confi_window ; l <= confi_window - 1; l ++ ){
		Posb.second = Pos2.second + l ;
		if ( sr_map[Posa].find(Posb) != sr_map[Posa].end() ){
		  vcf.sr += sr_map[Posa][Posb];
		}
	      }
	    }
	  }
	}
	else{ // SND
	  for (int k = - confi_window ; k < confi_window - 1; k ++ ){ 
	    Posa.second = Pos1.second + k ;
	    if ( un_map.find(Posa) != un_map.end() ){ 
	      vcf.sr += un_map[Posa];
	      vcf.pe += un_map[Posa];
	    }
	  }
	}
      }
    }

    for ( auto & i : container_vec[con].vcf_map ){
      for ( auto & vcf : i.second ){
	vcf.sr2 = vcf.sr;
      }
    }
  }

  std::cout << "Complete: SR\n";


  /////////////////////////////////////////////////////
  // Calculatiing PE
  for ( std::size_t con = 0 ; con < container_vec.size() ; con ++ ){
    for ( auto & i : container_vec[con].vcf_map ){
      Pos1 = i.first;
      Position Posa = Pos1;
      Posa.second += read_length / 2;
      Posa.second /= pe_window_size;
      for ( auto & vcf : i.second ){
	Pos2.first = id_ref_map[vcf.chr2];
	Pos2.second = vcf.pos2;
	Position Posb = Pos2;
	Posb.second += read_length / 2;
	Posb.second /= pe_window_size ;
	for (int k = - 2 * insert_size / pe_window_size ; k <= 2 * insert_size / pe_window_size ; k ++ ){
	  Posa.second = ( Pos1.second / pe_window_size ) + k ;
	  if ( pe_map.find(Posa) != pe_map.end() ){
	    for ( int l = - 2 * insert_size / pe_window_size ; l <= 2 * insert_size / pe_window_size ; l ++ ){
	      Posb.second = ( Pos2.second / pe_window_size ) + l ;
	      if ( pe_map[Posa].find(Posb) != pe_map[Posa].end() ){
		vcf.pe += pe_map[Posa][Posb];
	      }
	    }
	  }
	}
      }
    }
    for ( auto & i : container_vec[con].vcf_map ){
      for ( auto & vcf : i.second ){
	vcf.pe /= 2;
      }
    }

  // If sv type is SND, it is assumed that its PE is same with SR.
    for ( auto & i : container_vec[con].vcf_map ){
      for ( auto & vcf : i.second ){
	if ( vcf.svtype == "SND" ){
	  vcf.pe = vcf.sr;
	}
      }
    }
    
    for ( auto & i : container_vec[con].vcf_map ){
      for ( auto & vcf : i.second ){
	vcf.pe2 = vcf.pe;
      }
    }
  }

  std::cout << "Complete: PE \n";


  /////////////////////////////////////////////////////
  // Calculatiing ENTROPY
  std::map<Position,double> entropy_map;
  for ( auto & i : sr_map ){
    if ( i.second.size() < 2 ){
      entropy_map[i.first] = 0 ;
    }
    else{
      double sum = 0, prob = 0, entropy = 0 ;
      for ( auto & j : i.second ){
        sum += (double) j.second ;
      }
      for ( auto & j : i.second ){
        if ( sum > 0 ) prob = ( (double) j.second ) / sum ;
	if ( prob > 0 ) entropy = - prob * log ( prob );
      }
      if ( entropy < 0.000001 ) entropy = 0;
      entropy_map[i.first]=entropy;
    }
  }


  for ( std::size_t con = 0 ; con < container_vec.size() ; con ++ ){
    for ( auto & i : container_vec[con].vcf_map ){
      Pos1 = i.first;
      Position Posa = Pos1;
      for ( auto & vcf : i.second ){
	for (int k = - confi_window ; k < confi_window - 1; k ++ ){
	  Posa.second = Pos1.second + k ;
	  if ( entropy_map.find(Posa) != entropy_map.end() ){
	    vcf.entropy += entropy_map[Posa];
	  }
	}
	if ( vcf.svtype != "SND" ){
	  Pos2.first = id_ref_map[vcf.chr2];
	  Pos2.second = vcf.pos2;
	  Position Posb = Pos2;
	  for (int k = - confi_window ; k < confi_window - 1; k ++ ){
	    Posb.second = Pos2.second + k ;
	    if ( entropy_map.find(Posb) != entropy_map.end() ){
	      vcf.entropy2 += entropy_map[Posb];
	    }
	  }
	}
	else{
	  vcf.entropy2 = vcf.entropy;
	}
      }
    }
  }
  std::cout << "Complete: ENTROPY\n";



  /////////////////////////////////////////////////////
  // make info field
  for ( std::size_t con = 0 ; con < container_vec.size() ; con ++ ){
    for ( auto & i : container_vec[con].vcf_map ){
      for ( auto & j : i.second ){
	j.make_info(); // update info
      }
    }
  }
  std::cout << "Complete: feature calculation\n";

}
