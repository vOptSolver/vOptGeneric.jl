# vOptSolver: solver of Multiobjective linear optimization problems

The next release of vOptSolver (version 0.2) is scheduled for June 2017.
Currently the solver integrates exact algorithms for problems with two and three objectives.

[.] = forthcoming


Goals
- solver of multiobjective optimization problems for scientifics and practionners
- easy to formulate a problem, to provide data, to solve a problem, to collect the solutions, to analyze the outputs
- natural and intuitive for mathematicians, informaticians, engineers

Characteristics
- efficient and flexible
- free, open source, multi-platform, reusing existing specifications
- easy installation, no need of being expert in computer science

Purposes
- solving needs: methods and algorithms for performing numerical experiments
- research needs: support and primitives for the development of new algorithms
- pedagogic needs: environment for practicing of theories and algorithms

Background
- Julia programming language
- JuMP algebraic language
- Usual free and commercial MILP solvers (GLPK, CBC, CPLEX, GUROBI, etc.)

Problems / definition
- non-structured problems / algebraic language: ILP, MILP, LP
- structured problems / Application Programming Interface (API): OSP, LAP, [UKP, MKP, UDFLP, SSCFLP, UMFLP, CFLP, PATHS]

Algorithms
- the solving algorithms included compute an exact complete set of non-dominated ``points''
- generic algorithm for structured or non-structured discrete problems: Aneja-Nair method/2ILP, Chalmet's method/2ILP, e-constraint method/2ILP 
- specific algorithm for non-structured problem: branch-and-cut/2MILP
- specific algorithm for structured (MOCO/MOMILP) problem: 2OSP1980(jl), 2LAP2008 (C), [2UKP2009(C++), 2UDFLP2012 (C++), 2UMFLP2016(C++)]
- algorithms and datastructures for handling non-dominated points: AVL2016/2ILP(jl)

Inputs
- non-structured problems: direct with the languages (julia, jump), standard MOP format (ILP, MILP, LP), specific problem format (MILP)
- structured problems: direct with the languages (julia), specific problem format (2LAP, 2UKP, 2UFLP)

Outputs
- non-structured problems: standard 2MOP format (ILP, MILP, LP)
- structured problems: specific problem format (2LAP, 2UKP, 2UFLP)

Environment:
- free, open source under [GPL] licence
- julia 
- available; macOSX (12), linux (UBUNTU), windows (10)

Terms and acronyms used
- LP: Linear Program
- MILP: Mixed Integer Linear Program
- ILP: Integer Linear program
- MOCO: MultiObjective Combinatorial Optimization
- MOMILP: MultiObjective Mixed Integer Linear Program
- MOP: MultiObjective Program

- OSP: One machine Scheduling Problem
- LAP: Linear Assignment Problem
- UKP: Unidimensional Knapsack Problem
- MKP: Multidimensional Knapsack Problem
- UFLP: Uncapacitated Facility Location Problem
- UDFLP: Discrete Uncapacitated  Facility Location Problem
- SSCFLP: Single Source Capacitated Facility Location Problem
- UMFLP:  Mixed variables Uncapacitated Facility Location Problem
- CFLP: Capacitated Facility Location Problem
- PATHS: shortest paths problem

- Julia: name of the programming language
- JuMP: `Julia for Mathematical Optimization' is a modeling language for mathematical optimization embedded in Julia
- AVL tree is a self-balancing binary search tree
- API: `Application Programming Interface'
- GPL: `GNU General Public License'
- GLPK: `GNU Linear Programming Kit', an open source solver
- CBC: an open source solver
- CPLEX: a commercial solver
- GUROBI: a commercial solver


