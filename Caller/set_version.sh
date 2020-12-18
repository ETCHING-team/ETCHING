#!/usr/bin/env bash

version=$(cat ../version)
sed "s/etching_version=/etching_version=\"$version\"/g" my_vcf.cpp.0 > my_vcf.cpp

sed "s/const std::string version ;/const std::string version = \"$version\";/g" etching_caller.hpp.0 > etching_caller.hpp
