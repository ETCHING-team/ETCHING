#!/usr/bin/bash

if [ $# -ne 1 ]
then
    echo "Usage: checkerq.sh all/default"
    exit -1
fi

if [ $1 != "default" ]
then
    if [ $1 != "all" ]
    then
	echo "Usage: checkerq.sh all/default"
	exit -1
    fi
fi

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

if [ $1 == "all" ]
then
    # Check if cmake >=3.14 installed
    CMAKE_CHECK=$(cmake --version | head -n1 | awk '{print $3}' | awk -F "." '{if($1>=3 && $2>=14) print "YES"}')
    if [ $CMAKE_CHECK == "YES" ]
    then
	echo "make: OK"
    else
	echo "ERROR!!! Required: cmake >=3.14 for make all"
	exit -1
    fi
else
    # Check if cmake >=3.0 installed
    CMAKE_CHECK=$(cmake --version | head -n1 | awk '{print $3}' | awk -F "." '{if($1>=3) print "YES"}')
    if [ $CMAKE_CHECK == "YES" ]
    then
	echo "make: OK"
    else
	echo "ERROR!!! Required: cmake >=3"
	exit -1
    fi
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

if [ $1 == "default" ]
then
    if (( version < 6 ))
    then
	echo "Required: g++ >=6 for default"
	exit -1
    else
	echo -e "gcc/g++ >=6 for default: OK"    
    fi
else
    if [ $1 == "all" ]
    then
	if (( version < 7 ))
	then
	    echo "Required: g++ >=7 for full-installation"
	    exit -1
	else
	    echo -e "gcc/g++ >=7 for full-installation: OK"    
	fi
    else
	echo "ERROR!!! in checkreq.sh"
	exit -1
    fi
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
