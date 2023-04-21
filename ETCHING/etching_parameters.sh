#--------------------------------------------------------------------
# Copyright 2020- Bioinformatic and Genomics Lab.
# Hanyang University, Seoul, Korea
# Coded by Jang-il Sohn (sohnjangil@gmail.com)
#--------------------------------------------------------------------  

REQUIRED_SAMTOOLS_VERSION=1.13
PRESENT_DIR=$PWD
if [[ ${#ETCHING_PATH} == 0 ]]
then
    ETCHING_PATH=$(dirname $0)
    cd ${ETCHING_PATH}
    ETCHING_PATH=$PWD 2> /dev/null
    cd - 2> /dev/null
fi

ETCHING_INFO=
INCLUDE_PATH=${ETCHING_PATH}:${ETCHING_PATH}/../include:${CPATH}:${C_INCLUDE_PATH}:${CPLUS_INCLUDE_PATH}
if [[ ${#INCLUDE_PATH} != 0 ]]
then
    for SEARCH_PATH in $(echo ${INCLUDE_PATH} | sed "s/:/\n/g")
    do
	if [[ -f ${SEARCH_PATH}/etching_info.hpp ]]
	then
	    ETCHING_INFO=${SEARCH_PATH}/etching_info.hpp
	    break
	fi
    done
fi

if [[ ${#ETCHING_INFO} == 0 ]]
then
    echo -e "\033[31m## ------------------------------------"
    echo -e "## ERROR!!! We cannot find etching_info.hpp"
    echo -e "## ------------------------------------\033[0m"
    exit 1
fi

ETCHING_VERSION=$(grep "#define ETCHING_VERSION" ${ETCHING_INFO} | awk '{print $3}' | sed "s/\"//g")
RELEASE_DATE=$(grep "#define RELEASE_DATE" ${ETCHING_INFO} | awk -F "#define RELEASE_DATE " '{print $2}' | sed "s/\"//g")
MODEL_VERSION=$(grep "#define MODEL_VERSION" ${ETCHING_INFO} | awk '{print $3}' | sed "s/\"//g")
CONTACT=$(grep "#define ETCHING_CONTACT" ${ETCHING_INFO} | sed "s/#define ETCHING_CONTACT //" | sed "s/\"//g")

## Parameters
## required parameters
FIRST=
SECOND=
BAM=
CRAM=
GENOME=

PREFIX="etching"

THREADS=8

KMER_SIZE=31

WORKDIR_MES="same with -o"
WORKDIR=""

OUTDIR="./"

## gtf file for fusion-gene detection
## etching_fg_identifier requires this
ANNOTATION=

## k-mer parameters
KMER_CUTOFF=5

MAXK=10000

## sequencing read parameters
INSERTSIZE=500

## read orientation
ORIENT="FR"

## normal sequencing data as control
FIRST_CONT=
SECOND_CONT=
BAM_CONT=
CRAM_CONT=

## filter
FILTER=


### todo Update with --no-filter
## set 1 to skip etching_filter
## It will be 1, if no control (-1c, -2c, or -bc) and no k-mer set (-f or --kmer-table).
NO_FILTER=0

## pre-made k-mer table
KMER_TABLE=

## parameters for machine learning
ALGOR=
ALGOR_R=0
ALGOR_X=0

## default etching_sorter score
CUTOFF=0.4

## /path/to/ ETCHING machine learning models
ETCHING_ML_PATH=

## alignment tools
BWA=bwa
SAMTOOLS=samtools

## set KEEP_KMC=1 to keep KMC files for reuse
## else KEEP_KMC=0 to save storage
KEEP_KMC=0

## considering direction in etching_fg_identifier
STRAND_AWARE=0

## fusion-gene detecting window
FUSION_WINDOW=10000

## bwa mem option -T
BWA_T=30


## "default" for high recall
## "fast" for high speed
BAM_MODE="default"


## 
PREPROCESSING=1

SET_ENVIRONMENT_CHECK=0
CHECK_REQUIRED_CHECK=0
_LD_CHECK=0
