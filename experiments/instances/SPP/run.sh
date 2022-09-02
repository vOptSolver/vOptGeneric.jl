#!/bin/bash

for subfolder in ./BOSPP/*; do
    for file in $subfolder/*; do
        echo "$file"
        julia vOptBOSPP.jl "$file"
    done
done