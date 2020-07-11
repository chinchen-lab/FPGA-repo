#! /bin/bash

make clean
make

killall fpga

for i in $(seq 1 9)
do
    (../bin/./fpga $i > ../log/SRT_RR_0711_$i.log)&
done