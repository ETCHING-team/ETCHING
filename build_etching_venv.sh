#!/usr/bin/bash

# set environment
ETCHING_HOME=$PWD
ETCHING_BIN=${ETCHING_HOME}/bin
ETCHING_LIB=${ETCHING_HOME}/lib

# build etching_venv
python3 -m venv ${ETCHING_HOME}
laststatus=$?

if (( laststatus !=0 ))
then
    echo "ERROR!!! \"python3 -m venv etching_venv\" exited annormally."
    exit -1
fi

mv ${ETCHING_BIN}/activate ${ETCHING_BIN}/etching_venv
mv ${ETCHING_BIN}/activate.csh ${ETCHING_BIN}/etching_venv.csh
mv ${ETCHING_BIN}/activate.fish ${ETCHING_BIN}/etching_venv.fish

sed -i "s/VIRTUAL_ENV=.*/VIRTUAL_ENV=${ETCHING_HOME}/g" ${ETCHING_BIN}/etching_venv
sed -i "s/setenv VIRTUAL_ENV.*/setenv VIRTUAL_ENV ${ETCHING_HOME}/g" ${ETCHING_BIN}/etching_venv.csh
sed -i "s/set -gx VIRTUAL_ENV.*/set -gx VIRTUAL_ENV ${ETCHING_HOME} \$PATH/g" ${ETCHING_BIN}/etching_venv.fish

echo "export LD_LIBRARY_PATH=${ETCHING_LIB}:\$LD_LIBRARY_PATH" >> ${ETCHING_BIN}/etching_venv
echo "setenv LD_LIBRARY_PATH ${ETCHING_LIB}:\$LD_LIBRARY_PATH" >> ${ETCHING_BIN}/etching_venv.csh
echo "set -gx LD_LIBRARY_PATH ${ETCHING_LIB} \$LD_LIBRARY_PATH" >> ${ETCHING_BIN}/etching_venv.csh

# Install required Python modules with requirements.txt:
# pandas, numpy, scikit-learn (0.23.2), skranger (<=0.3), and xgboost
source ${ETCHING_BIN}/etching_venv
wget https://bootstrap.pypa.io/get-pip.py -O get-pip.py
python3 get-pip.py
pip3 install -r requirements.txt
pip3 list # Check if all requirements were properly installed
deactivate

exit 0
