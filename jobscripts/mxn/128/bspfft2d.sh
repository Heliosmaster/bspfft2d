#! /bin/bash
cd ../../..
make mxn
./fft2d 1 $1 $2 $3 $4
