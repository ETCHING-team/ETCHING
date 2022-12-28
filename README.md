# ETCHING-v1.4.1

### Efficient Detection of Chromosomal Rearrangements Using a Scalable k-mer Database of Multiple Reference Genomes and Variations

For raw FASTQ whole genome sequencing (WGS) data, preprocessing and alignment steps normally take about 12 hours (or more), 
and then SV prediction also take more than an hour (sometime days). ETCHING takes 2~3 hours for raw FASTQ (30X normal
and 50X tumor) on 30 threads on DELL 930 server.

Matched normal sequencing data was essentially required to predict somatic SVs from tumor sequencing data. 
However, ETCHING can predict somatic SVs without matched normal using a pan-genome k-mer set
(PGK2, http://big.hanyang.ac.kr/ETCHING/download.html) based on 894 human genomes and assemblies.

You can also find and download our codes, PGK2, and DEMO files from our website.

http://big.hanyang.ac.kr/ETCHING/

## Announcement

* Our manuscript was published in Nature Biomedical Engineering (https://doi.org/10.1038/s41551-022-00980-5).
* The next version (ETCHING-v1.4.2) is coming soon.
  * Comparing with the present version, it is 20% faster for FASTQ input.
  * Since ETCHING was developed focusing on FASTQ input, it was not fast for aligned BAM. Now, it can predict SVs in an hour. 
  * Additionally, we updated ETCHING to support CRAM. 

## Recent changes

#### v1.4.1
* Minor update (debug ETCHING/etching line 534)

#### v1.4.0

* Update with a graph theory of breakend, BND=(BP1,BP2), where BP=(chr,pos,dir).
* Machine learning modules for ranger and xgboost were replaced with cpp.
* Requirement changed
  * g++ >=6 (For full install, >=7)
  * cmake >=3 (For full install, >=3.14)

* Quality control module (fastp) was included in ```etching_filter```.
* New Pan-Genome K-mer set, PGK2 (based on 894 human genomes), was released (http://big.hanyang.ac.kr/ETCHING/download.html).

*See CHANGE.md for older updates.*




# Requirements

### System

* 64-bit LINUX with >=32GB RAM (at least >=16GB).
* Tested on Fedora workstation, Centos, and Ubuntu

### Software


* Required to compile
  * Default install: gcc, g++ (>=6.1.0), make, cmake (>=3.0), wget
  * Full install: gcc, g++ (>=7.1.0), make, cmake (>=3.14), wget

* Required to run
  * BWA, samtools



# Guide to ETCHING

We prepared a simple guide for CentOS/Fedora or Ubuntu/Debian/Mint users. You can skip this this step if all requirements were installed.

*Note: We tested this guide on Fedora32/33/34, CentOS7/8, Ubuntu16.04/18.04/20.04, Mint19/20, Debian11, and MX linux.*



## 1. Requirements

- #### CentOS/Fedora (or other Red Hat-based linux distros)

```bash
# Required programs 
sudo yum install -y gcc gcc-c++ make cmake bwa samtools wget
```

- #### Ubuntu/Debian/Mint (or other Debian-based distros)

```bash
## Required programs 
sudo apt install -y gcc g++ make cmake bwa samtools wget
```



## 2. Installation

Once, requirements were solved, you can install ETCHING as follows.

```bash
# Download ETCHING
git clone --depth=1 https://github.com/ETCHING-team/ETCHING.git

# Move to /path/to/ETCHING
cd ETCHING

# Default install (without XGBoost)
make

# Full install (g++ >=7 and cmake >=3.14 required)
make all

# Set your environment for ETCHING
echo "export ETCHING_HOME=$PWD" >> ~/.bashrc
echo "export PATH=$PWD/bin:\$PATH" >> ~/.bashrc
exec $SHELL
```




## 3. DEMO

```bash
# Change directory
cd /wherever/you/want/

# Download and decompress DEMO
wget http://big.hanyang.ac.kr/ETCHING/tiny_demo.tar.gz
# or you can find a demo file from our website
# http://big.hanyang.ac.kr/ETCHING/

# Decompress
tar zxvf tiny_demo.tar.gz
cd tiny_demo

# Run demo
etching -1 tumor_1.fq -2 tumor_2.fq -1c normal_1.fq -2c normal_2.fq -g Chr22.fa -f demo_PGK
```



# Pan-Genome k-mer set

The Pan-Genome k-mer(PGK) set is used to build PGK filter. Since we updated the PGK to PGK2, you can also download the PGK2 as below. 
If you have no matched-normal data, you are highly recommended to use the PGK2 in stead of PGK to call somatic SVs. 

```bash
# Move to etching directory
cd /somewhere/you/want/

# Download
wget http://big.hanyang.ac.kr/ETCHING/PGK2.tar.gz

# Decompress
tar zxvf PGK2.tar.gz
cd PGK2
ls PGK2

# Then, you must see PGK2.kmc_pre and PGK2.kmc_suf in PGK2
# Here, PGK2 is the name of k-mer set for ETCHING.
```

Alternatively, you can make your custom k-mer set.

In the absence of matched-normal sample, we recommand to use >=200 genomes to build your custom k-mer set exlcuding rare k-mers (<1% allele frequency). 
If you have 600 genomes, use ```-m 6``` option, which is the minimum frequency of k-mer to be included in your k-mer set.

```bash
make_pgk -i 600_genome.list -m 6 -o my_pgk
```
Here, ```600_genome.list``` is a file of fasta files of reference genomes.



# ETCHING on a ship (docker)

### Requirement

docker 



### Download docker image

```bash
# Download ETCHING docker image
wget http://big.hanyang.ac.kr/ETCHING/ETCHING_v1.4.0.docker.tar

# Load the image
docker load -i ETCHING_v1.4.0.docker.tar

# Check the image
docker images
```

You can see like this

|REPOSITORY|TAG|IMAGE ID|CREATED|SIZE|
|:---|:---|:---|:---|:---|
|etching|v1.4.0|63ffc48504f0|40 hours ago|3.26GB|



### Demo for docker user

Download our DEMO
```bash
# Download and decompress DEMO
wget http://big.hanyang.ac.kr/ETCHING/DEMO.tar.gz
tar zxvf DEMO.tar.gz
```

Run ETCHING with docker
```bash
docker run -i -t --rm -v /path/to/DEMO/:/work/ etching:1.4.0 etching \
-1 tumor_1.fq -2 tumor_2.fq -1c normal_1.fq -2c normal_2.fq -g small_genome.fa \
-a small_genome.gtf -f /work/demo_PGK -t 8
```
Here, ```etching:1.4.0``` is ```REPOSITORY``` and ```TAG``` of ETCHING docker image.

Replace ```/path/to/DEMO``` with ```/your/data/path/```.

Note: Keep ```/work/``` in the above command line.


Alternatively, you can run ETCHING inside docker container
```bash
docker run -i -t --rm -v /path/to/DEMO/:/work/ etching:1.4.0 /bin/bash

etching -1 tumor_1.fq -2 tumor_2.fq -1c normal_1.fq -2c normal_2.fq -g small_genome.fa \
-a small_genome.gtf -f /work/demo_PGK
```




# ETCHING on Amazon Web Service

ETCHING is also available on Amazon Web Service (AMI ID: ami-07c7a7d8934784df9; Region: us-east-1 (Northern Virginia)).
```bash
# Lunch EC2 instance from ETCHING AMI
# And connect to EC2 instance
ssh -i Your_Key ubuntu@Your_Instance_Address

# Decompress
tar zxvf ~/resources/DEMO.tar.gz 

# Run demo
cd ~/resources/DEMO
bash example.sh
```



# Related programs

### Filtration tool for PacBio long-reads
https://github.com/ETCHING-team/LR_Filter



### Benchmarking tool
https://github.com/ETCHING-team/etching_bench




# Contributors

Jang-il Sohn, Min-Hak Choi, Dohun Yi, A. Vipin Menon, and Jin-Wu Nam

Bioinformatics and Genomics Lab., Dept. of Life Science, Hanyang University, Seoul 04763, Korea




# Contact

#### Principal investigator

Jin-Wu Nam ([jwnam@hanyang.ac.kr]())



##### If he is unavailable, please email one of

Jang-il Sohn (sohnjangil@gmail.com)

Min-hak Choi (choiminhak1004@gmail.com)

Dohun Yi (kutarballoon@gmail.com)

----------------------------------------------------------------------------------
