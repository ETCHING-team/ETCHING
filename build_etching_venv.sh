#!/usr/bin/bash

# build etching_venv
ETCHING_HOME=$PWD
python3 -m venv ${ETCHING_HOME}
laststatus=$?

if (( laststatus !=0 ))
then
    echo "ERROR!!! \"python3 -m venv etching_venv\" exited annormally."
    exit -1
fi

mv bin/activate bin/etching_venv
mv bin/activate.csh bin/etching_venv.csh
mv bin/activate.fish bin/etching_venv.fish

sed -i "s/VIRTUAL_ENV=.*/VIRTUAL_ENV=\${ETCHING_HOME}/g" ${ETCHING_HOME}/bin/etching_venv
sed -i "s/setenv VIRTUAL_ENV.*/setenv VIRTUAL_ENV \${ETCHING_HOME}/g" ${ETCHING_HOME}/bin/etching_venv.csh
sed -i "s/set -gx VIRTUAL_ENV.*/set -gx VIRTUAL_ENV \${ETCHING_HOME} \$PATH/g" ${ETCHING_HOME}/bin/etching_venv.fish

echo "export LD_LIBRARY_PATH=${ETCHING_HOME}/lib:\$LD_LIBRARY_PATH" >> ${ETCHING_HOME}/bin/etching_venv
echo "setenv LD_LIBRARY_PATH ${ETCHING_HOME}/lib:\$LD_LIBRARY_PATH" >> ${ETCHING_HOME}/bin/etching_venv.csh
echo "set -gx LD_LIBRARY_PATH ${ETCHING_HOME}/lib \$LD_LIBRARY_PATH" >> ${ETCHING_HOME}/bin/etching_venv.csh

# Install required Python modules with requirements.txt:
# pandas, numpy, scikit-learn (0.23.2), skranger (<=0.3), and xgboost
source bin/etching_venv
wget https://bootstrap.pypa.io/get-pip.py -O get-pip.py
python3 get-pip.py
pip3 install -r requirements.txt
pip3 list # Check if all requirements were properly installed
deactivate

exit 0
