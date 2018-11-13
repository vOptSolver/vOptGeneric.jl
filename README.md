# vOptGeneric: part of vOptSolver for non-structured problems

[![Build Status](https://travis-ci.org/vOptSolver/vOptGeneric.jl.svg?branch=master)](https://travis-ci.org/vOptSolver/vOptGeneric.jl)
[![codecov.io](http://codecov.io/github/vOptSolver/vOptGeneric.jl/coverage.svg?branch=master)](http://codecov.io/github/vOptSolver/vOptGeneric.jl?branch=master)

**vOptSolver** is a solver of multiobjective linear optimization problems (MOCO, MOIP, MOMILP, MOLP).
This repository concerns **vOptGeneric**, the part of vOptSolver devoted to **multiobjective non-structured problems** (currently available: 2-IP, p-ILP). With vOptGeneric, the problem is expressed using JuMP algebraic language extended to multiple objectives. vOptGeneric runs on macOS, linux-ubuntu, and windows (local use), also on JuliaBox (distant use).

We suppose you are familiar with vOptSolver; if not, read first this [presentation](https://voptsolver.github.io/vOptSolver/).


## Instructions 
For a local use, a working version of:
- Julia must be ready; instructions for the installation are available [here](https://julialang.org/downloads/)
- your favorite MILP solver must be ready (GLPK is suggested); 
  instructions for the installation are available [here](http://www.juliaopt.org/JuMP.jl/0.18/installation.html)
  
### Run Julia

On linux or in the cloud (JuliaBox):

- open a console on your computer or in the cloud
- when the prompt is ready, type in the console `julia`

On macOS:

- locate the application `julia` and 
- click on the icon, the julia console comes to the screen

### Installation Instructions

Before your first local or distant use, 
1. run Julia and when the terminal is ready with the prompt `julia` on screen, 
2. add as follow the mandatory packages to your Julia distribution: 

```
julia> using Pkg
julia> Pkg.add("vOptGeneric")
julia> Pkg.add("GLPK") ; Pkg.add("GLPKMathProgInterface")
```

That's all folk; at this point, vOptGeneric is properly installed.

### Usage Instructions

When vOptGeneric is properly installed,

1. run Julia and when the terminal is ready with the prompt `julia` on screen, 
2. invoke vOptGeneric and the MILP solver to activate in typing in the console:
```
julia> using vOptGeneric
julia> using GLPK ; using GLPKMathProgInterface
```
vOptGeneric is ready. See examples for further informations and have fun with the solver! 

## Problems available

| Problem | Description                          | Output    | Method                       | Parameter (if required)  | Name          |
|:--------|:-------------------------------------|:---------:| ---------------------------: | ------------| :--------|
| 2-ILP   | bi-objective Integer Linear Program  | Y_N     | **:epsilon**                 | step = *realValue*       | ϵ-constraint  | 
| 2-ILP   | bi-objective Integer Linear Program  | Y_N     | **:chalmet** or **:Chalmet** | step = *realValue*       | Chalmet       |
| 2-ILP   | bi-objective Integer Linear Program  | Y_{SN}  | **:dicho** or **:dichotomy** | (none)                   | Aneja & Nair  |
| p-ILP | multi-objective Integer Linear Program | Y_{lex} | **:lex** or **:lexico**      | (none)                   | Lexicographic |


## Examples
The folder `examples` provides (1) source code of problems ready to be solved and (2) selected datafiles into different formats.

## Limitations
No special limitation; the solving strength of vOptGeneric is currently provided by the MILP solver (GLPK, Clp/Cbc, CPLEX, GUROBI, etc.) invoked.
