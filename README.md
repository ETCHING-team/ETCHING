---
# ETCHING

### Version 1.1.2a

### Efficient Detection of Chromosomal Rearrangements Using a Scalable k-mer Database of Multiple Reference Genomes and Variations

ETCHING takes 2-3 hours for WGS data with 30X normal and 50X tumor with 20 threads on DELL 930 server.
You can also find codes, k-mer set, and DEMO files in our website.

http://big.hanyang.ac.kr/ETCHING/

The demo is complete within 10 min on a desktop (AMD Ryzen 7 3700X 8-Core Processor).

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
	* [Docker](#docker)
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
git clone https://github.com/ETCHING-team/ETCHING.git
cd etching
make
echo "export PATH=$PWD/bin:$PATH" >> ~/.bashrc
echo "export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH" >> ~/.bashrc
echo "export ETCHING_ML_PATH=$PWD/ETCHING_ML_model" >> ~/.bashrc
```

or from our web server

```
wget http://big.hanyang.ac.kr/ETCHING/ETCHING_v1.1.2.tar.gz
tar zxvf ETCHING_v1.1.2.tar.gz
cd ETCHING_v1.1.2
make
echo "export PATH=$PWD/bin:$PATH" >> ~/.bashrc
echo "export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH" >> ~/.bashrc
echo "export ETCHING_ML_PATH=$PWD/ETCHING_ML_model" >> ~/.bashrc
```

### From executable file

Alternatively, you can just use our pre-compiled version of ETCHING.

```
wget https://github.com/ETCHING-team/ETCHING/releases/download/v1.1.2/ETCHING_v1.1.2_binary.tar.gz
tar zxvf ETCHING_v1.1.2_binary.tar.gz
cd ETCHING_v1.1.2_binary
echo "export PATH=$PWD/bin:$PATH" >> ~/.bashrc
echo "export ETCHING_ML_PATH=$PWD/ETCHING_ML_model" >> ~/.bashrc
```

or

```
wget http://big.hanyang.ac.kr/ETCHING/ETCHING_v1.1.2_binary.tar.gz
tar zxvf ETCHING_v1.1.2_binary.tar.gz
cd ETCHING_v1.1.2_binary
echo "export PATH=$PWD/bin:$PATH" >> ~/.bashrc
echo "export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH" >> ~/.bashrc
echo "export ETCHING_ML_PATH=$PWD/ETCHING_ML_model" >> ~/.bashrc
```

### From docker image

Download the docker image from our web (http://big.hanyang.ac.kr/ETCHING/download.html)
And load ETCHING docker image
```
wget http://big.hanyang.ac.kr/ETCHING/etching_docker_v1.1.2.tar
docker load -i etching_v1.1.2.tar
```

Check if the docker image is loaded properly

```
docker images
```

Output should be like below

|REPOSITORY|TAG|IMAGE ID|CREATED|SIZE|
|:---|:---|:---|:---|:---|
|etching|1.1.2|16647cac9a99|40 hours ago|3.5GB|


### Pan-Genome k-mer set

After installation of ETCHING, you need pan-genome k-mer set (PGK).

```
cd /path/to/etching
wget http://big.hanyang.ac.kr/ETCHING/PGK.tar.gz
tar zxvf PGK.tar.gz
```

Then, you will see the files ```PGK_20200103.kmc_pre``` and ```PGK_20200103.kmc_suf``` 
in the directory ```PGK```:

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

Keep the order and ```-k31```, ```-ci1```, and ```-fa``` options. The ```genome_list.txt``` is a file of
reference genomes' list, which must be separated line-by-line. Use any name you prefer for ```your_kmer_db```. 
If you use the same same, you will see ```your_kmer_db.pre``` and ```your_kmer_db.suf``` in the directory.

## Usage

After installation, you can download and run demo

```
wget http://big.hanyang.ac.kr/ETCHING/DEMO.tar.gz
tar zxvf DEMO.tar.gz
cd DEMO
etching -i demo.conf
```

If you want the list of options, check with this command

```
etching -h
```


### Example of a config file

You can write a config file as follows, if you want to predict SVs from tumor samples (tumor_1/2) with 
matched normal samples (normal_a/b_1/2), and ready-made k-mer database (PGK).

```
#####################################################################
#                                                                           
# sample     : (tumor) sample in fastq                                      
# filter_fq  : (matched normal) control sample as filter in fastq 
# filter_fa  : Assembly, genome, or normal samples as filter in fasta
# filter_db  : Pre-made k-mer set (prefix of database), such as PGK
# genome     : Reference genome for mapping 
# mapper     : bwa or minimap2 
# annotation : for fusion-gene identifying in gtf format 
# 
# NOTE!!!       
# Keep the type names.              
# Use filter_db only once, if you needed it.
#####################################################################


# sample (or tumor) sequencing files
# ETCHING supports gzip (fq.gz) files, 
# but it will slow down 3-5 times in the Filter module.
sample  tumor_1.fastq
sample  tumor_2.fastq

# filter (or normal) samples
filter_fq       normal_a_1.fastq
filter_fq       normal_a_2.fastq

# Additional fastq files can be added
filter_fq       /path/to/normal_b_1.fastq
filter_fq       /path/to/normal_b_2.fastq

# Ready-made k-mer database made with KMC3
filter_db       /path/to/PGK_20200103

# If you want to add more k-mers from 
# assemblies (or genomes) to k-mer set, 
# use these options
filter_fa       Ash1.7.fa
filter_fa       ../../T2T.fasta
filter_fa       /path/to/CHM1.fa

# BWA (or Minimap2) indexed reference genome
genome          /path/to/hg19.fa

# Set one of these if you need to specify mapper.
# default is bwa
#mapper         /path/to/bwa
#mapper         /path/to/minimap2

# samtools
# set this if you need to specify it
#samtools       /path/to/samtools

# An annotation file in gtf format for FG-identifier
# Use this option to predic FGs.
# Note: only gtf is supported yet.
annotation      hg19.gtf
```

If you have no matched normal data, PGK must be helpful to select tumor specific reads.

### Example execution 

#### WGS data

If you want to run ETCHING with 30 threads (```-t 30```) for sequencing data of length 150 bp 
(```-L 150```), run the following command:

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

### Docker

In case of using ETCHING docker image,

```
docker run -i -t --rm -v /path/to/DEMO/:/work/ etching:1.1.2 etching -i demo.conf
```
Mount your directory containing all required files, such as PGK, reference fasta, sample and normal 
fastq  to '/work/' directory in the docker container. The mount point on docker '/work/' should not be
changed. Note that file paths in a configure file should point the mounted files in the docker container. 
So path should be either relative path starting from your mount point or absolute path starting from 
'/work/' directory. To say, if you have the data path '/path/to/mounted/directory/sample.fq', and mounted
'/path/to/mounted/directory' to '/work/' on docker container, then the file path on configure file should
be either 'sample.fq' or '/work/sample.fq'.

Alternatively, you can run ETCHING inside docker container
```
docker run -i -t --rm -v /local/path/to/example/directory/:/work/ etching:1.1.2 /bin/bash
etching -i example.conf
```
The ETCHING binary and ML models are located on '/etching/'



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
