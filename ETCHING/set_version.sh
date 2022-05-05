#!/usr/bin/env bash

version_line=$(grep "ETCHING" ../version)
ml_version=$(grep "ML_version" ../version | awk '{print $2}')

version=$(grep "ETCHING" ../version | awk '{print $1}')
release=$(grep "ETCHING" ../version | awk '{print $2}')

awk -F "=" -v version=$version -v release=$release '{if($1=="VERSION") print "VERSION=\""version" "release"\"" ; else print $0}' etching | awk -F "=" -v ml_version=$ml_version '{if($1=="MODEL_VERSION") print "MODEL_VERSION="ml_version ; else print $0}' > tmp
mv tmp etching
