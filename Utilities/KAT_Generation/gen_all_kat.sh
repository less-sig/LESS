#!/usr/bin/env bash

### uncomment the following lines if you also want to build the project
#mkdir -p build 
#cd build 
#cmake ..
#make 
#cd ..

for i in build/*
do
    if [[ -x "${i}" ]]
    then
        echo Generating KATs for ${i}
        ./${i}
    fi
done

mkdir -p KAT
mv *.req ./KAT
mv *.rsp ./KAT
