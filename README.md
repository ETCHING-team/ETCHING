---
# ETCHING

### Version 1.1.5

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
### From docker image

Download the docker image from our web (http://big.hanyang.ac.kr/ETCHING/download.html)
And load ETCHING docker image
```
wget http://big.hanyang.ac.kr/ETCHING/etching_docker_v1.1.3.tar
docker load -i etching_v1.1.3.tar
```

Check if the docker image is loaded properly

```
docker images
```

Output should be like below

|REPOSITORY|TAG|IMAGE ID|CREATED|SIZE|
|:---|:---|:---|:---|:---|
|etching|1.1.3|16647cac9a99|40 hours ago|3.5GB|


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
If you have no matched normal data, PGK must be helpful to select tumor specific reads.
./ㅐㅔ;.,

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
etching -1 tumor_data/tumor1.fq -2 tumor_data/tumor2.fq \
-1c normal_data/normal1.fastq -2c normal_data/normal2.fastq \
-g small_genome/small_genome.fa \
-a small_genome/small_genome.gtf \
-p OUTPUT \
-f kmer_set/demo \
-m /path/to/ETCHING_ML_model
```

If you want the list of options, check with this command

```
etching -h
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

If you use `-a annotation.gtf` option, ETCHING will
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
etching_filter -1 tumor_1.fq -2 tumor_2.fq \
-1c normal_1.fastq -2c normal_2.fastq \
-p OUTPUT \
-t 30 \
-g genome.fa \
-p OUTPUT \
-f PGK \
-m /path/to/ETCHING_ML_model

```

For panel, us `-T P` option

```
etching_filter -1 tumor_1.fq -2 tumor_2.fq \
-1c normal_1.fastq -2c normal_2.fastq \
-p OUTPUT \
-t 30 \
-g genome.fa \
-p OUTPUT \
-f PGK \
-m /path/to/ETCHING_ML_model \
 P -A
```


2. Caller

```
etching_caller -b filtered_read.sort.bam -g genome.fa -o OUTPUT
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
docker run -i -t --rm -v /path/to/DEMO/:/work/ etching:1.1.3 etching -i demo.conf
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
docker run -i -t --rm -v /local/path/to/example/directory/:/work/ etching:1.1.3 /bin/bash
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
