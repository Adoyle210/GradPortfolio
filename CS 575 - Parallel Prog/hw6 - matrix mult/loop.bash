#!/bin/bash

for s in 32 64 128 256 512 1024
do
        for b in 1 2 4 8 16
        do
                g++ -DMATW=$s -DLOCALSIZE=$b -o proj06 proj06.cpp /usr/local/apps/cuda/cuda-10.1/lib64/libOpenCL.so.1.1 -lm -fopenmp
                ./proj06
        done
done