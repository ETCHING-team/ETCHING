//--------------------------------------------------------------------
// Copyright 2020. Bioinformatic and Genomics Lab.
// Hanyang University, Seoul, Korea
// Coded by Jang-il Sohn (sohnjangil@gmail.com)
//--------------------------------------------------------------------

#include "find_fg.hpp"

const int FF = 0;
const int FR = 1;
const int RF = 2;
const int RR = 3;

void read_bp_pair(const std::string BP_pair_file, std::map < int , std::map < Position , std::set < Position > > > & BP_pair)
{
  std::ifstream fin(BP_pair_file.c_str());

  int smap[4];
  smap[0]=3;
  smap[1]=1;
  smap[2]=2;
  smap[3]=0;
  
  std::string chr1;
  std::string chr2;
  int pos1;
  int pos2;
  std::string Str;
  int strand;

  Position Pos1;
  Position Pos2;

  while ( fin >> chr1 >> pos1 >> chr2 >> pos2 >> Str ){
    if ( Str == "FF" ) strand = FF;
    else if ( Str == "FR" ) strand = FR;
    else if ( Str == "RF" ) strand = RF;
    else if ( Str == "RR" ) strand = RR;
    Pos1=std::make_pair(chr1,pos1);
    Pos2=std::make_pair(chr2,pos2);
    BP_pair[strand][Pos1].insert(Pos2);
    strand = smap[strand];
    BP_pair[strand][Pos2].insert(Pos1);
  }
  fin.close();
}


void read_gene_list(const std::string gene_list_file,
		    std::map < std::string , std::map < Loci , std::string > > & gene_list_for,
		    std::map < std::string , std::map < Loci , std::string > > & gene_list_rev)
{
  std::string id;
  std::size_t start;
  std::size_t end;
  std::string strand;
  std::string tmp;

  std::ifstream fin (gene_list_file.c_str());
  while ( fin >> id >> start >> end >> strand >> tmp ){
    if ( strand == "+" ){
      gene_list_for[id][std::make_pair(start,end)]=tmp;
    }
    else{
      gene_list_rev[id][std::make_pair(start,end)]=tmp;
    }
  }

  fin.close();
}



void find_fusion_gene ( std::map < int , std::map < Position , std::set < Position > > > & BP_pair,
			std::map < std::string , std::map < Loci , std::string > > & gene_list_for,
			std::map < std::string , std::map < Loci , std::string > > & gene_list_rev,
			std::set < std::string > & fusion_genes)
{
  std::string chr1;
  std::string chr2;
  int pos1;
  int pos2;

  std::string output;

  Position bp1;
  Position bp2;

  // FF
  for ( auto & i : BP_pair[FF] ){
    for ( auto & j : i.second ){
      bp1=i.first;
      bp2=j;
      std::set < std::string > gene1;
      std::set < std::string > gene2;
      chr1=bp1.first ; pos1=bp1.second ;
      chr2=bp2.first ; pos2=bp2.second ;
      for ( auto & j : gene_list_for[chr1] ){
	if ( pos1 >= j.first.first && pos1 <= j.first.second ){
	  gene1.insert(j.second);
	}
      }
      for ( auto & j : gene_list_for[chr2] ){
	if ( pos2 >= j.first.first && pos2 <= j.first.second ){
	  gene2.insert(j.second);
	}
      }
      std::string g1;
      std::string g2;

      std::set < std::string > gene1_uniq;
      std::set < std::string > gene2_uniq;
      for ( auto & g : gene1 ) if ( gene2.find(g) == gene2.end() ) gene1_uniq.insert(g);
      for ( auto & g : gene2 ) if ( gene1.find(g) == gene1.end() ) gene2_uniq.insert(g);
      gene1 = gene1_uniq;
      gene2 = gene2_uniq;

      for ( auto & g : gene1 ){
	g1+=g+"/";
      }
      for ( auto & g : gene2 ){
	g2+=g+"/";
      }

      if ( g2.size() == 0 ){
	g2 += "*/";
      }

      if ( g1.size() > 0 && g2.size() > 0 && g1 != g2){
	g1.pop_back();
	g2.pop_back();
	output = g1 + "\t" + g2 + "\t" + chr1 + "\t" + std::to_string(pos1) + "\t+\t" + chr2 + "\t" + std::to_string(pos2) + "\t+";
	fusion_genes.insert(output);
      }
    }
  }


  // FR
  for ( auto & i : BP_pair[FR] ){
    for ( auto & j : i.second ){
      bp1=i.first;
      bp2=j;
      std::set < std::string > gene1;
      std::set < std::string > gene2;
      
      chr1=bp1.first ; pos1=bp1.second;
      chr2=bp2.first ; pos2=bp2.second;
      
      for ( auto & j : gene_list_for[chr1] ){
	if ( pos1 >= j.first.first && pos1 <= j.first.second ){
	  gene1.insert(j.second);
	}
      }
      for ( auto & j : gene_list_rev[chr2] ){
	if ( pos2 >= j.first.first && pos2 <= j.first.second ){
	  gene2.insert(j.second);
	}
      }
    
      std::string g1;
      std::string g2;

      std::set < std::string > gene1_uniq;
      std::set < std::string > gene2_uniq;
      for ( auto & g : gene1 ) if ( gene2.find(g) == gene2.end() ) gene1_uniq.insert(g);
      for ( auto & g : gene2 ) if ( gene1.find(g) == gene1.end() ) gene2_uniq.insert(g);
      gene1 = gene1_uniq;
      gene2 = gene2_uniq;

      for ( auto & g : gene1 ){
	g1+=g+"/";
      }
      for ( auto & g : gene2 ){
	g2+=g+"/";
      }

      if ( g2.size() == 0 ){
	g2 += "*/";
      }

      if ( g1.size() > 0 && g2.size() > 0 && g1 != g2 ){
	g1.pop_back();
	g2.pop_back();
	output = g1 + "\t" + g2 + "\t" + chr1 + "\t" + std::to_string(pos1) + "\t+\t" + chr2 + "\t" + std::to_string(pos2) + "\t-";
	fusion_genes.insert(output);
      }
    }
  }

  // RF
  for ( auto & i : BP_pair[RF] ){
    for ( auto & j : i.second ){
      bp1=i.first;
      bp2=j;
      
      std::set < std::string > gene1;
      std::set < std::string > gene2;
      
      chr1=bp1.first ; pos1=bp1.second;
      chr2=bp2.first ; pos2=bp2.second;
      
      for ( auto & j : gene_list_rev[chr1] ){
	if ( pos1 >= j.first.first && pos1 <= j.first.second ){
	  gene1.insert(j.second);
	}
      }
      for ( auto & j : gene_list_for[chr2] ){
	if ( pos2 >= j.first.first && pos2 <= j.first.second ){
	  gene2.insert(j.second);
	}
      }
      
      std::string g1;
      std::string g2;
      
      std::set < std::string > gene1_uniq;
      std::set < std::string > gene2_uniq;
      for ( auto & g : gene1 ) if ( gene2.find(g) == gene2.end() ) gene1_uniq.insert(g);
      for ( auto & g : gene2 ) if ( gene1.find(g) == gene1.end() ) gene2_uniq.insert(g);
      gene1 = gene1_uniq;
      gene2 = gene2_uniq;
      
      for ( auto & g : gene1 ){
	g1+=g+"/";
      }
      for ( auto & g : gene2 ){
	g2+=g+"/";
      }
      
      if ( g2.size() == 0 ){
	g2 += "*/";
      }
      
      if ( g1.size() > 0 && g2.size() > 0 && g1 != g2 ){
	g1.pop_back();
	g2.pop_back();
	output = g1 + "\t" + g2 + "\t" + chr1 + "\t" + std::to_string(pos1) + "\t-\t" + chr2 + "\t" + std::to_string(pos2) + "\t+";
	fusion_genes.insert(output);
      }
    }
  }
  // RR
  for ( auto & i : BP_pair[RR] ){
    for ( auto & j : i.second ){
      bp1=i.first;
      bp2=j;

      std::set < std::string > gene1;
      std::set < std::string > gene2;

      chr1=bp1.first ; pos1=bp1.second;
      chr2=bp2.first ; pos2=bp2.second;

      for ( auto & j : gene_list_rev[chr1] ){
	if ( pos1 >= j.first.first && pos1 <= j.first.second ){
	  gene1.insert(j.second);
	}
      }
      for ( auto & j : gene_list_rev[chr2] ){
	if ( pos2 >= j.first.first && pos2 <= j.first.second ){
	  gene2.insert(j.second);
	}
      }

      std::string g1;
      std::string g2;

      std::set < std::string > gene1_uniq;
      std::set < std::string > gene2_uniq;
      for ( auto & g : gene1 ) if ( gene2.find(g) == gene2.end() ) gene1_uniq.insert(g);
      for ( auto & g : gene2 ) if ( gene1.find(g) == gene1.end() ) gene2_uniq.insert(g);
      gene1 = gene1_uniq;
      gene2 = gene2_uniq;

      for ( auto & g : gene1 ){
	g1+=g+"/";
      }
      for ( auto & g : gene2 ){
	g2+=g+"/";
      }

      if ( g2.size() == 0 ){
	g2 += "*/";
      }

      if ( g1.size() > 0 && g2.size() > 0 && g1 != g2 ){
	g1.pop_back();
	g2.pop_back();
	output = g1 + "\t" + g2 + "\t" + chr1 + "\t" + std::to_string(pos1) + "\t-\t" + chr2 + "\t" + std::to_string(pos2) + "\t-";
	fusion_genes.insert(output);
      }
    }
  }
}


int find_fg ( std::string BP_pair_file, std::string gene_list_file){
  std::map < std::string , std::map < Loci , std::string > > gene_list_for;
  std::map < std::string , std::map < Loci , std::string > > gene_list_rev;

  std::map < int , std::map < Position , std::set < Position > > > BP_pair;

  std::set < std::string > fusion_genes;

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Reading SV vcf file
  //

  read_bp_pair(BP_pair_file,BP_pair);


  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Reading gene list
  //

  read_gene_list(gene_list_file,gene_list_for,gene_list_rev);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Find fusion-genes
  //

  find_fusion_gene(BP_pair,gene_list_for,gene_list_rev,fusion_genes);

  for ( auto & i : fusion_genes ){
    std::cout << i << "\n";
  }

  return 0;
}
