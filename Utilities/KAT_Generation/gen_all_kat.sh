#!/bin/bash

# mkdir -p build 
# cd build 
# cmake ..
# make 
# cd ..

for i in build/bin/*
do
    echo Generating KATs for $i
    ./$i
done

mv *.req ../../KAT/
mv *.rsp ../../KAT/
