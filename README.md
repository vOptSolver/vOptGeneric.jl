# **vOptSolver**: solver of Multiobjective linear optimization problems

vOptSolver is currently supported by the ANR/DFG-14-CE35-0034-01 research project [(link)](https://voptproject.wordpress.com/). 

The next release of vOptSolver (version 0.2) is scheduled for June 2017.
It integrates exact algorithms for computing a complete set of non-dominated points for structured and non-structured optimization problems with two [and three] objectives ([ ] = forthcoming).

Coordinator: Prof. Dr. Xavier Gandibleux [(contact)](http://www.univ-nantes.fr/gandibleux-x)

## Presentation

### Aims
- Solver of multiobjective linear optimization problems for scientifics and practionners
- Easy to formulate a problem, to provide data, to solve a problem, to collect the outputs, to analyze the solutions
- Natural and intuitive use for mathematicians, informaticians, engineers

### Purposes
- Solving needs: methods and algorithms for performing numerical experiments
- Research needs: support and primitives for the development of new algorithms
- Pedagogic needs: environment for practicing of theories and algorithms

### Characteristics
- Efficient, flexible, evolutive solver
- Free, open source, multi-platform, reusing existing specifications
- Easy installation, no need of being expert in computer science

### Background
- Julia programming language [(link)](http://julialang.org/)
- JuMP algebraic language [(link)](https://jump.readthedocs.io/en/latest/)
- Usual free and commercial MILP solvers (GLPK, CPLEX, GUROBI)

## Features

### Problems / Definition
- Non-structured problems / algebraic language: 
    -  ILP: Integer Linear Program
    -  MILP: Mixed Integer Linear Program
    -  LP: Linear Program
- Structured problems / Application Programming Interface (API): 
    -  OSP: One machine Scheduling Problem
    -  LAP: Linear Assignment Problem 
    -  [UKP, MKP, UDFLP, SSCFLP, UMFLP, CFLP, PATHS]

### Algorithms
The solving algorithms included compute an exact complete set of non-dominated points
- Generic algorithm for structured or non-structured discrete problems: 
    - Aneja & Nair method / 2ILP
    - Chalmet's et al. method / 2ILP
    - epsilon-constraint method / 2ILP 
- Specific algorithm for non-structured problem: 
    - branch-and-cut / 2MILP
- Specific algorithm for structured (MOCO/MOMILP) problem: 
    - 2OSP1980 (Julia)
    - 2LAP2008 (C)
    - [2UKP2009 (C++), 2UDFLP2012 (C++), 2UMFLP2016 (C++)]
- Algorithms and datastructures for handling non-dominated points: 
    - AVL2016/2ILP (Julia)

### Inputs
- Non-structured problems: 
    - direct with the provided languages (Julia, JuMP)
    - standard MOP format (ILP, MILP, LP)
    - specific problem format (MILP)
- Structured problems: 
    -  direct with the languages (Julia), 
    -  specific problem format (2LAP, 2UKP, 2UFLP)

### Outputs
- Non-structured problems: 
    - standard 2MOP format (ILP, MILP, LP)
- Structured problems: 
    - specific problem format (2LAP, 2UKP, 2UFLP)

## Running

### Information
- Julia is available on macOS, linux, windows for a local use or online on [JuliaBox](https://juliabox.com/) for a distant use
- vOptSolver is free, open source under [GPL] licence, tested with Julia 0.5 on macOS and linux-Ubuntu

### Installation
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

### Example
It is time to try vOptSolver:

```
julia>  
```

---

Terms and acronyms used
- LP: Linear Program
- MILP: Mixed Integer Linear Program
- ILP: Integer Linear program
- MOCO: MultiObjective Combinatorial Optimization
- MOMILP: MultiObjective Mixed Integer Linear Program
- MOP: MultiObjective Program
- OSP: One machine Scheduling Problem
- LAP: Linear Assignment Problem
- UKP: Unidimensional 01 Knapsack Problem
- MKP: Multidimensional 01 Knapsack Problem
- UFLP: Uncapacitated Facility Location Problem
- UDFLP: Discrete Uncapacitated  Facility Location Problem
- SSCFLP: Single Source Capacitated Facility Location Problem
- UMFLP:  Mixed variables Uncapacitated Facility Location Problem
- CFLP: Capacitated Facility Location Problem
- PATHS: shortest paths problem
- Julia: name of the programming language
- JuMP: stands for Julia for Mathematical Optimization, a modeling language for mathematical optimization embedded in Julia
- AVL tree is a self-balancing binary search tree
- API: stands for Application Programming Interface
- GPL: stands for GNU General Public License
- GLPK: stands for GNU Linear Programming Kit, an open source solver
- CPLEX: a commercial solver
- GUROBI: a commercial solver


