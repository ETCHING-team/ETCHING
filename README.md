# ETCHING

### Version 1.4.0

### Efficient Detection of Chromosomal Rearrangements Using a Scalable k-mer Database of Multiple Reference Genomes and Variations

ETCHING takes about 3 hours for WGS data with 30X normal and 50X tumor on 30 threads on DELL 930 server.
You can also find codes, k-mer set, and DEMO files in our website.

http://big.hanyang.ac.kr/ETCHING/

 

## Recent changes

#### v1.4.0

* Update with a graph theory of breakend, BND=(BP1,BP2), where BP=(chr,pos,dir).
* Machine learning modules for ranger and xgboost were replaced with cpp.
* Requirement changed
  * g++ >=6 (For full install, >=7)
  * cmake >=3 (For full install, >=3.14)

* Quality control module (fastq) was included in etching_filter.
* New Pan-Genome K-mer set, PGK2 (based on 894 human genomes), was released (http://big.hanyang.ac.kr/ETCHING/).

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

If you have no matched-normal data, you are highly recommended to use the pan-genome k-mer set2 (PGK2) to call somatic SVs. 

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

Alternatively, you can make your own k-mer set.

We recommand >=200 genomes to remove rare k-mers (<1% AF). 
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
wget http://big.hanyang.ac.kr/ETCHING/ETCHING_v1.4.0.docker.saved.tar

# Load the image
docker load -i ETCHING_v1.4.0.docker.saved.tar

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
docker run -i -t --rm -v /path/to/DEMO/:/work/ etching:1.4.0 etching -1 tumor_1.fq -2 tumor_2.fq -1c normal_1.fq -2c normal_2.fq -g small_genome.fa -a small_genome.gtf -f /work/demo_PGK -t 8
```
Here, ```etching:1.4.0``` is ```REPOSITORY``` and ```TAG``` of ETCHING docker image.

Replace ```/path/to/DEMO``` with ```/your/data/path/```.

Note: Keep ```/work/``` in the above command line.


Alternatively, you can run ETCHING inside docker container
```bash
docker run -i -t --rm -v /path/to/DEMO/:/work/ etching:1.4.0 /bin/bash

etching -1 tumor_1.fq -2 tumor_2.fq -1c normal_1.fq -2c normal_2.fq -g small_genome.fa -a small_genome.gtf -f /work/demo_PGK
```


----------------------------------------------------------------------------------


# Related programs

### Filtration tool for PacBio long-reads
https://github.com/ETCHING-team/LR_Filter

### Benchmarking tool
https://github.com/ETCHING-team/etching_bench


----------------------------------------------------------------------------------


# Contributors

Jang-il Sohn, Min-Hak Choi, Dohun Yi, A. Vipin Menon, and Jin-Wu Nam

Bioinformatic and Genomics Lab., Dept. of Life Science, Hanyang University, Seoul 04763, Korea

----------------------------------------------------------------------------------


# Contact

If you have any issues, please contact us

   Jang-il Sohn (sohnjangil@gmail.com)

   Jin-Wu Nam (jwnam@hanyang.ac.kr)

----------------------------------------------------------------------------------
