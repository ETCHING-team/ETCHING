---
# ETCHING

### Version 1.1.1

### Efficient Detection of Chromosomal Rearrangements Using a Scalable k-mer Database of Multiple Reference Genomes and Variations

ETCHING takes 2-3 hours for WGS data with 30X normal and 50X tumor with 20 threads.

---

## Table of contents

  * [Requirement](#requirement)
  * [Installation](#installation)
    * [From source code](#from-source-code)
    * [From executable file](#from-executable-file)
    * [Pan-Genome k-mer set](#pan-genome-k-mer-set)
  * [Usage](#usage)
    * [Example of a config file](#example-of-a-config-file)
    * [Pan-genome k-mer - PGK](#pan-genome-k-mer---pgk)
    * [Example execution](#example-execution)
    * [Step-by-step execution](#step-by-step-execution)
  * [Contributors](#contributors)
  * [Contact](#contact)


## Requirement

### System

* 64-bit LINUX with >=64GB RAM (at least >=16GB).

* C++11 with g++>=4.7.0

* python3 with pandas, numpy, scikit-learn, skranger, and xgboost packages. You can simply install the packages as follows:
```
pip3 install pandas numpy scikit-learn skranger xgboost
```
* BWA (https://github.com/lh3/bwa.git) or alternatively Minimap2 (https://github.com/lh3/minimap2.git)
 
* samtools (http://www.htslib.org/)

You can install ETCHING even if you did not install the python packages, but ETCHING-Sorter is not working properly.

##### We tested this version on CentOS 7.4.1708 with g++ >4.7 and python-3.6.10.


## Installation

### From source code
Download code from git
```
git clone https://github.com/sohnjangil/etching.git
cd etching
make
echo "export PATH=$PWD/bin:$PATH" >> ~/.bashrc
echo "export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH" >> ~/.bashrc
echo "export ETCHING_ML_PATH=$PWD/ETCHING_ML_model" >> ~/.bashrc
```
or from our web server
```
wget http://big.hanyang.ac.kr/ETCHING/ETCHING_v1.1.1.tar.gz
tar zxvf ETCHING_v1.1.1.tar.gz
cd ETCHING_v1.1.1
make
echo "export PATH=$PWD/bin:$PATH" >> ~/.bashrc
echo "export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH" >> ~/.bashrc
echo "export ETCHING_ML_PATH=$PWD/ETCHING_ML_model" >> ~/.bashrc
```

### From executable file
Alternatively, you can just use our pre-compiled version of ETCHING.
```
wget https://github.com/ETCHING-team/ETCHING/releases/download/v1.1.1/ETCHING_v1.1.1_binary.tar.gz
tar zxvf ETCHING_v1.1.1_binary.tar.gz
cd ETCHING_v1.1.1_binary
echo "export PATH=$PWD/bin:$PATH" >> ~/.bashrc
echo "export ETCHING_ML_PATH=$PWD/ETCHING_ML_model" >> ~/.bashrc
```
or
```
wget http://big.hanyang.ac.kr/ETCHING/ETCHING_v1.1.1_binary.tar.gz
tar zxvf ETCHING_v1.1.1_binary.tar.gz
cd ETCHING_v1.1.1_binary
echo "export PATH=$PWD/bin:$PATH" >> ~/.bashrc
echo "export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH" >> ~/.bashrc
echo "export ETCHING_ML_PATH=$PWD/ETCHING_ML_model" >> ~/.bashrc
```



### Pan-Genome k-mer set

After installation of ETCHING, you may need pan-genome k-mer set (PGK) especially when you have no matched normal of tumor samples.
```
cd /path/to/etching
wget http://big.hanyang.ac.kr/ETCHING/PGK.tar.gz
tar zxvf PGK.tar.gz
```
Then, you will see the files ```PGK_20200103.kmc_pre``` and ```PGK_20200103.kmc_pref``` in the directory ```PGK```:
```
ls PGK/
PGK_20200103.kmc_pre
PGK_20200103.kmc_suf
```
Here, ```PGK_20200103``` is the name of k-mer set to be used for ETCHING.

Alternatively, you can make your own k-mer set as follows:
```
kmc -k31 -v -ci1 -fa @genome_list.txt your_kmer_db ./
```
Keep the order and ```-k31```, ```-ci1```, and ```-fa``` options. The ```genome_list.txt``` is a file of reference genomes' list, which must be separated line-by-line. Use any name you prefer for ```your_kmer_db```. If you use the same same, you will see ```your_kmer_db.pre``` and ```your_kmer_db.suf``` in the directory.

## Usage

You can run ETCHING with a config file
```
etching -i exmple.conf
```

### Example of a config file

You can write a config file as follows, if you want to predict SVs from tumor samples (tumor_1/2) with matched normal samples
(normal_a/b_1/2), reference genomes (hg19, T2T, and KOREF), and ready-made k-mer database (PGK).
```
# sample (or tumor) sequencing files
sample  tumor_1.fastq
sample  tumor_2.fastq

# filter (or normal) samples
filter_fq       normal_a_1.fastq
filter_fq       normal_a_2.fastq

filter_fq       /path/to/normal_b_1.fastq
filter_fq       /path/to/normal_b_2.fastq

# genome (in fasta format) to be used as filter
filter_fa       hg19.fasta
filter_fa       ../../T2T.fasta
filter_fa       /path/to/KOREF.fasta

# Ready-made k-mer database made with KMC3
filter_db       /path/to/PGK_20200103

# Indexed reference genome
genome          hg19.fa

# set this if you need to specify it
# default is bwa
#mapper         /path/to/bwa
#mapper         /path/to/minimap2

# samtools
# set this if you need to specify it
#samtools       /path/to/samtools

# An annotation file in gtf format for FG-identifier
# Use this option to predic FGs.
annotation      hg19.gtf
```

If you have no matched normal data, PGK must be helpful to select tumor specific reads.

### Example execution 

#### WGS data

If you want to run ETCHING with 30 threads (```-t 30```) for sequencing data of length 150 bp (```-L 150```), run the following command:
```
etching -i example.conf -t 30 -L 150 -p OUTPUT
```

#### Panel data

In the case of targeted panel data (```-T P```), we recommend to use all split-reads with ```-A``` option:
```
etching -i example.conf -t 30 -L 150 -T P -A -p OUTPUT
```

#### Output

The above command will give you two vcf files:
```
OUTPUT.BND.etching_sorter.vcf
OUTPUT.SV.etching_sorter.vcf
```

In OUTPUT.BND.etching_sorter.vcf, all vcf lines are recoded in breakend (BND) level.
In OUTPUT.SV.etching_sorter.vcf, you will see SV typs, such as DEL, DUP, INV, or TRA. 

#### Fusion-gene

If you use a line of annotation file in config file, ETCHING will
make two more files of fusion-genes:
```
OUTPUT.BND.fusion_gene.txt
OUTPUT.SV.fusion_gene.txt
```

### Step-by-step execution

If you want to run ETCHING step-by-step, follow the steps.

1. Filter

For WGS,
```
etching_filter -i example.conf -t 30
```

For panel,
```
etching_filter -i example.conf -t 30 -T P -A
```

Note: The configure file (```example.conf```) and ```-T P``` option are available only here.

2. Caller
```
etching_caller -b filtered_read.sort.bam -g /path/to/hg19.fa -o OUTPUT
```
This command returns two output, ```OUTPUT.BND.vcf``` and ```OUTPUT.SV.vcf```.

3. Sorter
```
etching_sorter -i OUTPUT.BND.vcf -o OUTPUT.BND
etching_sorter -i OUTPUT.SV.vcf -o OUTPUT.SV
```
4. FG_identifiter
```
etching_fg_identifier output.BND.etching_sorter.vcf hg19.annot.gtf > output.BND.FG.txt
etching_fg_identifier output.SV.etching_sorter.vcf hg19.annot.gtf > output.SV.FG.txt
```

----------------------------------------------------------------------------------
## Contributors

Min-Hak Choi, Jang-il Sohn, Dohun Yi, A. Vipin Menon, and Jin-Wu Nam

Bioinformatic and Genomics Lab., Hanyang University, Seoul 04763, Korea

----------------------------------------------------------------------------------
## Contact

If you have any issues, please contact us freely.

   Jang-il Sohn (sohnjangil@gmail.com)

   Jin-Wu Nam (jwnam@hanyang.ac.kr)

----------------------------------------------------------------------------------