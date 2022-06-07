#!/bin/bash

for file in ./*.jl; do
    echo "$file"
    julia "$file"
done

# for file in ./*.jl; do
#     echo "$file"
#     julia "$file" 0.5
# done

# for file in ./*.jl; do
#     echo "$file"
#     julia "$file" 1
# done