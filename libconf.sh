version=$(g++ --version | head -n 1 | awk '{print $3}' | awk -F "." '{print $1}')

libversion=

if [ $version -lt 4 ]
then
    echo "[Requirement]"
    echo "You need gcc >4.7.0"
    exit -1
fi

if [ $version -eq 4 ]
then
    subversion=$(g++ --version | head -n 1 | awk '{print $3}' | awk -F "." '{print $2}')
    if [ $subversion -lt 7 ]
    then
	echo "[Requirement]"
	echo "You need gcc >4.7.0"
	exit -1
    else
	libversion="4.7.0"
    fi
fi


if [ $version -eq 5 ]
then
    libversion="5.1.0"
fi

if [ $version -eq 6 ]
then
    libversion="6.1.0"
fi

if [ $version -eq 7 ]
then
    libversion="7.1.0"
fi

if [ $version -eq 8 ]
then
    libversion="8.1.0"
fi

if [ $version -eq 9 ]
then
    libversion="9.1.0"
fi

if [ $version -eq 10 ]
then
    libversion="10.1.0"
fi

cd lib
ln -sf ${libversion}/*.so ./
cd ..
