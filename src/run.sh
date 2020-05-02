#! /bin/bash

make clean
make

#killall fpga

for i in $(seq 1 9)
do
    (../bin/./fpga $i > ../log/out_times_$i.log )&
done