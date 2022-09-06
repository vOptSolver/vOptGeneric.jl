#!/bin/bash

for subfolder in ./BOSPA/*.txt; do
    echo "$file"
    julia vOptBOSPA.jl "$file"
done

julia latexWriter.jl