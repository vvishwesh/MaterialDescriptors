#!/bin/bash

# Building MaterialDescriptors
if [ -d ../lib ]; then
    if [ -d lib ]; then
        cp -r ../lib/*.jar lib/
    else
        cp -r ../lib .        
    fi 
fi


find ../src/materialdescriptors/ -name *.java > javafiles.txt
javac -cp lib/commons-io-2.8.0.jar @javafiles.txt -encoding utf-8 -d .

if [ "$?" != "0" ]; then
    rm javafiles.txt
	echo "Failed to create MaterialDescriptors.jar."
    exit -1
fi

rm javafiles.txt


echo "Manifest-Version: 1.0" > manifest.mf
echo "Main-Class: materialdescriptors.MaterialDescriptors" >> manifest.mf
echo "Class-Path: lib/commons-io-2.8.0.jar" >> manifest.mf
echo >> manifest.mf

jar cvfm MaterialDescriptors.jar manifest.mf materialdescriptors 

if [ "$?" = "0" ]; then
    rm -rf manifest.mf materialdescriptors
else
	echo "Failed to create MaterialDescriptors.jar."
    exit -1
fi

echo "--------------------- Done building MaterialDescriptors.jar ---------------------"
