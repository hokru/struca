# structure analysis on the command line

This is a **simple** tool to help compare structures and compute rajectory data on the command line.


## Getting started

### How?

* edit FC in Makefile (default gfortran)
* run `struca -h` to see all options


#### available features:
* quaternion RMSD & alignment (all-atom)
* comparison of all unique internal coordinates  + statistik
* comparison of H-bonds 
* comparison of Rotational constants
* request custom primitives for analysis (bonds/angles/torsions)
* print internal coordinates, H-bonds, Rotational constants for single molecule
* [soon] trajectory analysis

### prerequisites

- **BLAS/LAPACK**: performance not critical, any library will suffice 
- **Fortran compiler**: developed and tested with gfortran (4.7).
- **git**: git is needed to include the git version in the output.

### scientific references
- RMSD quaternion fit: DOI:10.1002/jcc.20110
- cite code as (suggestion): H.Kruse www.github/hokru/struca V. *git version*
