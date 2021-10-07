#!/usr/bin/bash

# Check if awk is installed
check_awk=$(which awk)
if [ ${#check_awk} == 0 ]
then
    echo "Warning: awk is required to run ETCHING"
else
    echo "awk: OK"
fi


# Check if make is installed
check_make=$(which make)
if [ ${#check_make} == 0 ]
then
    echo "ERROR!!! Required: make"
    exit -1
else
    echo "make: OK"
fi

# check g++ and its version
check_gcc=$(which gcc)
if [ ${#check_gcc} == 0 ]
then
    echo "ERROR!!! Required: gcc"
    exit -1
else
    echo "gcc: OK"
fi

check_gpp=$(which g++)
if [ ${#check_gpp} == 0 ]
then
    echo "ERROR!!! Required: g++"
    exit -1
else
    echo "g++: OK"
fi


version=0
version=$(g++ --version | head -n 1 | awk '{print $3}' | awk -F "." '{print $1}')

if (( version < 4 ))
then
    echo "Required: g++ >4.7.0"
else
    if (( version == 4 ))
    then
	subversion=$(g++ --version | head -n 1 | awk '{print $3}' | awk -F "." '{print $2}')
	if (( subversion < 7 ))
	then
	    echo "ERROR!!! Required: g++ >4.7.0"
	    exit -1
	else
	    echo -e "gcc/g++: OK"
	fi
    else
	echo -e "gcc/g++: OK"
    fi
fi


# Check python3 version
pv=$(python3 --version | awk '{print $2}')
p3_1=$(echo $pv | awk -F "." '{print $1}')
p3_2=$(echo $pv | awk -F "." '{print $2}')


if (( p3_1 != 3 ))
then
    echo "ERROR!!! Required: Python3 (3.6, 3.7, or 3.8)"
    echo "Your python version: $pv"
    exit -1
fi

if (( p3_2 != 6 && p3_2 != 7 && p3_2 != 8))
then
    echo "ERROR!!! Required: Python3 (3.6, 3.7, or 3.8)"
    echo "Your python version: $pv"
    exit -1
else
    echo -e "Python3 (3.6-3.8): OK"
fi


# check python3-venv
message=$(python3 -m venv etching_venv)
rm -rf etching_venv
check_venv=$(echo $message | grep "Failing command:" | wc -l)
if (( check_venv !=0 ))
then
    echo "ERROR!!! \"python3 -m venv etching_venv\" exited abnormally. Please check it. "
    echo "----------------------"
    echo $message
    exit -1
else
    echo -e "python3-venv: OK"
fi


# Check if bwa is installed
check_bwa=$(which bwa)
if [ ${#check_bwa} == 0 ]
then
    echo "Warning: bwa is required to run ETCHING"
else
    echo "bwa: OK"
fi

# Check if samtools is installed
check_samtools=$(which samtools)
if [ ${#check_samtools} == 0 ]
then
    echo "Warning: samtools is required to run ETCHING"
else
    echo "samtools: OK"
fi


exit 0
