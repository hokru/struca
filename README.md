# structure analysis on the command line

This is a **simple** tool to help compare structures and compute trajectory data on the command line.


## Getting started

### How?

* edit FC in Makefile (default gfortran)
* run `struca -h` to see all options


# available options:
`    options:`

`    -h                       this help`

`    -comp <ref. structure> <structure-to-compare>    compare 2 molecules`

`    -struc <structure>                               analyse singluar molecule`


`    -bthr/-athr/-tthr        bond/angle/torsion thresholds`


`  <structure> formats: XMOL TMOL `


### prerequisites

- **BLAS/LAPACK**: performance not critical, any library will suffice 
- **Fortran compilter**: developed and tested with gfortran (4.7).
- **git**: git is needed to include the git version in the output.

### scientific references
- RMSD quaternion fit: DOI:10.1002/jcc.20110
- cite code as (suggestion): H.Kruse www.github/hokru/struca V. *git version*
