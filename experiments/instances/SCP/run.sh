#!/bin/bash

for file in ./BOSCP/*.dat; do
    echo "$file"
    julia vOptBOSCP.jl "$file"
done