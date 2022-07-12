#!/bin/bash

for file in ./*.jl; do
    echo "$file" 
    julia "$file"
done
