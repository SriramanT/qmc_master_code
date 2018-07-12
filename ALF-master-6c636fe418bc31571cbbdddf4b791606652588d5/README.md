# ALF (Algorithms for Lattice Fermions) #
## General information ##
This version of the **A**lgorithms for **L**attice **F**ermions package provides a general code for the finite temperature auxiliary field Quantum Monte Carlo algorithm.       The code  is engineered to  be able simulate any model that can be written in terms of  sums of single body operators, of squares of single body operators and single body operators coupled to an Ising field with  given dynamics. We  provide predefined types that allow  the user to specify the model, the  Bravais lattice  as well as equal time and time displaced observables.     The code supports an MPI implementation.   Examples such as the Hubbard model on the Honeycomb lattice  as well as the Hubbard model  on the square lattice coupled to a transverse Ising field are  provided and discussed in the [documentation](https://alf.physik.uni-wuerzburg.de/ALF/ALF/raw/master/Documentation/ALF_v1.0.pdf).  For questions and issues you can reach us at [alf@physik.uni-wuerzburg.de](mailto:alf@physik.uni-wuerzburg.de).

You can [Download ALF-v1.0.zip](https://alf.physik.uni-wuerzburg.de/ALF/ALF/repository/archive.zip?ref=master) or sign up to help us improve the package, create issues etc.

The Hamiltonians we can consider read:
![The Hamiltonian0](https://alf.physik.uni-wuerzburg.de/ALF/ALF/raw/master/Images/Hamiltonian0.png)
 where
![The Hamiltonian1](https://alf.physik.uni-wuerzburg.de/ALF/ALF/raw/master/Images/Hamiltonian1.png)

Here Z denotes an Ising spin variable with predefined dynamics. If your model can be written in this form then it will be amenable to the ALF. 
## PREREQUISITES ##
ALF is Fortran2003 conforming and should be compilable by every recent Fortran compiler with the help of a lapack library.

Libraries: LAPACK and BLAS. Regularly tested are MKL, OpenBLAS and the netlib.org implementation.

MPI support is optional.

Compiler: gfortran, ifort, PGI compiler

We have collected some hints for various unices [here](Installation.md).

## CONFIGURATION FOR COMPILATION ##
Basic configuration of this package happens directly in the top-level Makefile.
Here you can provide compiler information and flags as well as the location of LAPACK and Blas.
A simple 
```
. ./setenv.sh
make
```
should do the trick to compile the libraries, the analysis programs and the Example program.
## FILES AND DIRECTORIES ##

**Libraries**    Libraries.   

**Prog**   Main program and subroutines. The command **make** will generate the **Examples.out** executable. 

**Analysis** Analysis programs. The command **make** will generate three executables required  for the error analysis  of scalar, equal time and time displaced observables. 

**Start**   This directory contains the files required to start a run. In particular it contains the parameter file   that specifies the model the lattice and various   parameters for the Monte Carlo run and  error analysis. 

**Examples** This directory provides a set of short example runs.   

**Documentation**  We have included in the file [ALF_v1.0.pdf](https://alf.physik.uni-wuerzburg.de/ALF/ALF/raw/master/Documentation/ALF_v1.0.pdf)
an extensive documentation.

## TUTORIAL ##

A Tutorial can be found under https://alf.physik.uni-wuerzburg.de/ALF/ALF-Tutorial

## TESTING ##

We have about 30 tests that test various parts of the program in the folder testsuite.
As testing framework we employ CTest.
From the subfolder testsuite the tests can be run as follows
```
- mkdir tests
- cd tests
- cmake ..
- make
- make test
```

## CONTRIBUTORS ##

F. Assaad, M. Bercx, S. Beyl, F. Goth, J. Hofmann, M. Hohenadler, F. Parisen Toldin, T. Sato, J. Schwab and Z. Wang.

## CITATION ##

The documentation is available on the arxiv at [arXiv:1704.00131 [cond-mat.str-el]](https://arxiv.org/abs/1704.00131).

If you use the ALF package, please mention  the following in the acknowledments:

The auxillary field QMC simulations were carried out with the ALF package available at https://alf.physik.uni-wuerzburg.de .

## LICENSE ##
The various works that make up the ALF project are placed under licenses that put
a strong emphasis on the attribution of the original authors and the sharing of the contained knowledge.
To that end we have placed the ALF source code under the GPL version 3 license (see license.GPL and license.additional)
and took the liberty as per GPLv3 section 7 to include additional terms that deal with the attribution
of the original authors(see license.additional).
The documentation of the ALF project by the ALF contributors is licensed under a Creative Commons Attribution-ShareAlike 4.0 International License (see Documentation/license.CCBYSA)
We mention that we link against parts of lapack which licensed under a BSD license(see license.BSD).


[Disclaimer](https://www.uni-wuerzburg.de/sonstiges/impressum/) [Data Policy](dsgvo.md)