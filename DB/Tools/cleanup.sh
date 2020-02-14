#!/bin/bash

cd ..
cd DB/ReflectionGroups

for i in {5..5}
do	
    cd G"$i"_CHEVIE/Cherednik
    
    mkdir -p GGOR/Generic
    git mv Generic/GGOR/* GGOR/Generic
    
    cd ..
    cd ..
    
    #mkdir -p Cherednik/Generic/GGOR/CenterGenerators
    #regex='^CherednikGenericGGORCenterGenerator[0-9]+'
    #for f in CherednikGenericGGORCenterGenerator*; do
    #	test -f "$f" || continue
    #	if [[ $f =~ $regex ]]
    #	then
    #		g=$(echo ${f} | sed 's/CherednikGenericGGORCenterGenerator//')
    #		git mv $f Cherednik/Generic/GGOR/CenterGenerators/$g
    #	fi
    #done
    
done