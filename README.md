This repository accompagnies the Wiley paper "Generalization Riccati Recursion". The purpose of the paper is dual:

- reproducable benchmark results and possibility to test out performance on other platforms and for other problems
  
- possible re-use of code, specifically the gen_riccati part of this repo


This repository contains


- gen_riccati: a library that implements the generalization of Riccati recursion from the paper, the only dependency here is blasfeo

dependencies: blasfeo

- symspals: a library named that allows setting up sparse linear systems symbolically, the library is interfaced to three sparse linear solvers: mumps, ma57 and pardiso

dependencies: sparse linear solvers: mumps, ma57 and pardiso

- the code for performing the benchmark of the paper: here we have a dependency on fatrop for setting up the problems

dependencies: symspals, gen_riccati and fatrop


references:



Known issues: not possible to use ma57 and mumps with metis library at the same time because ma57 