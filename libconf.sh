version=$(g++ --version | head -n 1 | gawk '{print $3}' | gawk -F "." '{print $1}')

libversion=

if (( version < 4 ))
then
    echo "[Requirement))"
    echo "You need gcc >4.7.0"
    exit -1
fi

if (( version == 4 ))
then
    subversion=$(g++ --version | head -n 1 | gawk '{print $3}' | gawk -F "." '{print $2}')
    if (( $subversion < 7 ))
    then
	echo "[Requirement))"
	echo "You need gcc >4.7.0"
	exit -1
    else
	libversion="4.7.0"
    fi
fi

if (( version > 10 ))
then
    version=10
fi

libversion="${version}.1.0"

cd lib
ln -sf ${libversion}/*.so ./
cd ..
