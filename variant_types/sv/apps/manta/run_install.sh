#!/usr/bin/env bash

wget https://github.com/Illumina/manta/archive/refs/tags/v1.6.0.tar.gz
gunzip v1.6.0.tar.gz
tar -xf v1.6.0.tar
cd manta-1.6.0
mkdir build
cd build
../configure
make -j
