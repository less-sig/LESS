#!/usr/bin/env bash

VALGRIND=valgrind

mkdir -p build 
cd build 
cmake .. > /dev/null 2>&1
make -j > /dev/null 2>&1
cd ..

rm -r ./build/*.massif

for file in build/bin/*
do
	bfile="$(basename ${file})"
	echo ${bfile}
	${VALGRIND} --tool=massif --detailed-freq=1 --max-snapshots=1000 --stacks=yes --main-stacksize=100000000 --max-stackframe=184121216 --massif-out-file=build/${bfile}.massif ./${file} > /dev/null 2>&1
done
