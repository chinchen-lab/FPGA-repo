#! /bin/bash

make clean
make

killall fpga_subch

for i in $(seq 1 9)
do
    (../bin/./fpga_subch $i > ../log/subch$i.log )&
done