#!/bin/bash

for i in build/bin/*
do
    echo Generating KATs for $i
    ./$i
done

mv *.req ../../KAT/
mv *.rsp ../../KAT/
