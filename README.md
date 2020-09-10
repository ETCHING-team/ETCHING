---
# ETCHING

### Version 1.0.1

### Efficient Detection of Chromosomal Rearrangements Using a Scalable k-mer Database of Multiple Reference Genomes and Variations

ETCHING takes 2-3 hours for WGS data with 30X normal and 50X tumor.

### Note!!!

If you are using v1.0 (initial buggy version), please replace with the version >=1.0.1.

---
## Table of contents

  * [Requirement](#requirement)
  * [Installation](#installation)
  * [Usage](#usage)
    * [Example of a config file](#example-of-a-config-file)
    * [Pan-genome k-mer - PGK](#pan-genome-k-mer---pgk)
    * [Example execution](#example-execution)
    * [Step-by-step execution](#step-by-step-execution)
  * [Contributors](#contributors)
  * [Contact](#contact)
---

## Requirement

### System

* 64-bit LINUX with >=64GB RAM (at least >=16GB).

* C++11 with g++>=4.7.0

* python2.7 with pickle, pandas, numpy, xgboost, and scikit-learn packages.

* bwa mem (or Minimap2)

* samtools

(Python3 will be supported soon)

You can install ETCHING even if you did not install the python packages, but ETCHING-Sorter is not working properly.

##### We tested this version on CentOS 7.4.1708 with g++-4.8.5 and python 2.7.15.

----------------------------------------------------------------------------------
## Installation

Download, compile, and set path
```
git clone https://github.com/sohnjangil/etching.git
cd etching

make

export PATH=/path/to/etching/bin:$PATH
export LD_LIBRARY_PATH=/path/to/etching/lib:$LD_LIBRARY_PATH
export ETCHING_ML_PATH=/path/to/etching/ETCHING_ML_model
```
You can add this line in bashrc or bash_profile.
```
echo "export PATH=/path/to/etching/bin:$PATH" >> ~/.bashrc
echo "export LD_LIBRARY_PATH=/path/to/etching/lib:$LD_LIBRARY_PATH" >> ~/.bashrc
echo "export ETCHING_ML_PATH=/path/to/etching/ETCHING_ML_model" >>/.bashrc
```
If you want to install ETCHING for all users, copy all files as follows.
```
sudo cp /path/to/etching/bin/* /usr/bin
sudo cp /path/to/etching/lib/* /usr/lib
```
or to /usr/local/bin, and /usr/local/lib.

----------------------------------------------------------------------------------
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

# control (or normal) samples
filter_fq       normal_a_1.fastq
filter_fq       normal_a_2.fastq

filter_fq       /path/to/normal_b_1.fastq
filter_fq       /path/to/normal_b_2.fastq

filter_fa       hg19.fasta
filter_fa       ../../genome/T2T.fasta
filter_fa       /path/to/KOREF.fasta

# Ready-made k-mer database made with KMC3
filter_db       /path/to/somewhere/PGK

# Indexed reference genome
genome          hg19.fa

# set this if you need to specify it
# default is bwa
#mapper         /path/to/bwa
#mapper         /path/to/minimap2

# samtools
# set this if you need to specify it
#samtools       /path/to/samtools

# annotation file in gtf format for FG-identifier
# if you did not use this, ETCHING will not predic FGs.
#annotation      hg19.gtf
```

### Pan-genome k-mer - PGK

We prepared our PGK database on our website.
You can download or make your own k-mer set of reference genomes using KMC3 instead of our PGK
```
kmc -k31 -v -ci1 -fa @genome_list.txt your_kmer_db ./
```
or update
```
kmc_tools simple /path/to/PGK -ci1 union your_kmer_db -ci1 PGK_updated
```

### Example execution

Run this command, if you want to run ETCHING with the example.conf file sequencing reads of length 150 bp (-L 150) targeted panel
data (-T P) on 30 threads (-t 30) 
```
etching -i example.conf -t 30 -L 150 -T P -p OUTPUT
```

The above command will give you two vcf files:
```
OUTPUT.BND.etching_sorter.vcf
OUTPUT.SV.etching_sorter.vcf
```

In BND.vcf, all vcf lines are recoded in breakend (BND) level without SV type,
such as DEL, DUP, INV, or TRA. You can see SV typs in SV.vcf file.

If you add (uncomment) a line of annotation file in config file, it will
make two more files of fusion-genes:
```
OUTPUT.BND.fusion_gene.txt
OUTPUT.SV.fusion_gene.txt
```

### Step-by-step execution

If you want to run ETCHING step-by-step,

1. Filter
```
etching_filter -i example.conf -t 30 -T P 
```
2. Caller
```
etching_caller -b filtered_read.sort.bam -g /path/to/hg19.fa -o output
```
3. Sorter
```
etching_sorter -i output.BND.vcf -o output.BND
etching_sorter -i output.SV.vcf -o output.SV
```
4. FG_identifiter
```
etching_fg_identifier output.BND.etching_sorter.vcf hg19.annot.gtf > output.BND.fusion_gene.txt
etching_fg_identifier output.SV.etching_sorter.vcf hg19.annot.gtf > output.SV.fusion_gene.txt
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