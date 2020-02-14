#!/bin/bash

# Sets links for degrees of complex reflection groups

cd ..
cd GrpMat

for i in {4..37}
do
    #cd G"$i"_Magma_Dual
    #ln -s ../G"$i"_Magma/Degrees.m Degrees.m
    #cd ..
    cd G"$i"_Magma
    rm ../G"$i"_CHEVIE/Degrees.m
    cp Degrees.m ../G"$i"_CHEVIE
    rm ../G"$i"_CHEVIE_Dual/Degrees.m
    cp Degrees.m ../G"$i"_CHEVIE_Dual
    rm ../G"$i"_Magma_Dual/Degrees.m
    cp Degrees.m ../G"$i"_Magma_Dual
    #rm Representations.m
    #git add CharacterData.m
    cd ..
    #cd G"$i"_CHEVIE_Dual
    #ln -s ../G"$i"_Magma/Degrees.m Degrees.m
    #cd ..
    #cd G"$i"_CHEVIE
    #mv Representations.m Representations_0.m
done