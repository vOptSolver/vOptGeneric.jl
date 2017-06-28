# Solver of multiobjective linear optimization problems

This repository concerns **vOptGeneric**, the part of **vOptSolver** devoted to **multiobjective non-structured problems**. With vOptGeneric, the problem is expressed using JuMP algebraic language extended to multiple objectives.

We suppose you are familiar with vOptSolver; if not, read first this [presentation](https://voptsolver.github.io/vOptSolver/).


## Instructions 

### Information
- Julia is available on macOS, linux, windows for a local use or online on [JuliaBox](https://juliabox.com/) for a distant use
- vOptSolver is free, open source under [GPL] licence, tested with Julia 0.5 on macOS and linux-Ubuntu

### Installation Instructions
For a local use, a working version of:
- Julia must be ready; instructions for the installation are available [here](https://julialang.org/downloads/)
- your favorite MILP solver must be ready (mandatory for dealing with non-structured problems; GLPK is suggested); 
  instructions for the installation are available [here](http://jump.readthedocs.io/en/latest/installation.html)

Before your first local or distant use, run Julia and when the terminal is ready, add as follow the two mandatory packages to your Julia distribution: 


```
julia> Pkg.add("JuMP")
julia> Pkg.clone("vOptSolver")
```

That's all folk! 

### Usage Instructions

Install it, and then have fun with vOptSolver. 

Resources availables:

- [Examples](https://github.com/xgandibleux/vOptSolver/tree/master/examples/): 
    - problems ready to be solved
- [Documentation](https://github.com/xgandibleux/vOptSolver/tree/master/doc/): 
    - description of the Application Programming Interface 
