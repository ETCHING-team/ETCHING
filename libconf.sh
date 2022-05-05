version=$(g++ --version | head -n 1 | awk '{print $3}' | awk -F "." '{print $1}')

gcc_version=

if (( version < 4 ))
then
    echo "[Requirement))"
    echo "You need gcc >4.7.0"
    exit -1
fi

if (( version == 4 ))
then
    subversion=$(g++ --version | head -n 1 | awk '{print $3}' | awk -F "." '{print $2}')
    if (( subversion < 7 ))
    then
	echo "[Requirement]"
	echo "You need gcc >4.7.0"
	exit -1
    else
	gcc_version="4.7.0"
    fi
else
    if (( version > 10 ))
    then
	version=10
    fi
    gcc_version="${version}.1.0"
fi



cd lib

for i in filter caller fg_identifier
do
    if [ ! -f libetching_${i}.so ]
    then
	ln -sf ${gcc_version}/libetching_${i}.so ./ 
    fi
done

cd ..

exit 0
