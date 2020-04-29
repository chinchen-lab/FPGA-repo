#! /bin/bash
cd ../src
make clean
make
cd ../bin
./fpga 2
./fpga 1
