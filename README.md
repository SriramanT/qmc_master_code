# qmc_master_code

This repository is part of a MSc in Physics Engineering final project consisting of applying Quantum Monte Carlo (QMC) to a problem of interacting fermions. We refer the interested reader to the thesis available at

[Thesis](https://github.com/fmonteir/qmc_master_thesis/blob/master/thesis/thesis.pdf)

The aim of this work was to carry out a theoretical study (with particular emphasis on numerical aspects) of the properties of a TMD (transition metal dichalcogenide) nanoribbon - a graphene-like 2D nanostructure where electron interactions are particularly relevant - using a QMC method.

## Quantum Monte Carlo for Interacting Fermions

Capturing the effects of electron correlations is not an easy task. The difficulty lies in devising a numerical method to solve the many-body Schrödinger equation in a reasonable amount of computer time. The development and application of unbiased methods is a central point in correlated electron systems, particularly in 2D. Naive methods have exponential complexity in the system size, which motivates an approach based on Monte Carlo sampling.

Quantum Monte Carlo (QMC) methods are among the few unbiased methods available to date. In general, they circumvent the exponential complexity hurdle making it algebraic instead. However, for fermionic systems, a sign oscillation deems the algorithm exponential in the system size and inverse temperature, hence not very effective. This is due to the antisymmetric nature of the many-fermion wavefunction. One of the main tasks of this project was to investigate whether this issue would impede the simulation of a minimal model of 2D nanostructures made out of novel graphene-like 2D materials.

Some approximations exist to deal with the sign problem, and they perform differently depending on the problem at hand. However, it might be the case that the sign problem is either absent, as it happens for a class of models, or is not very serious, still allowing an unbiased and accurate simulation of the system at hand. Going beyond the fermion sign problem barrier in quantum simulations is an open topic of research and many possibilities have been put forward recently, namely an approach based on neural networks (see thesis and references therein).

In short, this masters thesis is about the implementation of an algorithm to deal with tight-binding problems for 2D interacting electronic models. Ultimately, the goal was to write a code that would simulate a specific interacting electron system: a transition metal dichalcogenide (TMD) nanoribbon.

TMD's are graphene-like 2D materials that are promising from both a theoretical and an application perspective. A nanoribbon is a 2D  nanostructure that is much longer on one direction than on the other (like a ribbon), so that electronic edge states become relevant and lead to unusual properties. In practice, this means that one can use periodic boundary conditions on the longer direction, and open boundary conditions on the other.

An example of an interesting property of these nanostructures that we investigate is magnetism. Furthermore, one might be interested in the different phases that arise within the system and in how do the transitions between them occur. For example, recent papers point at the possibility of topological superconductivity in TMD nanoribbons. This is a many-body effect that is only captured numerically by state of the art techniques such as QMC.

Our work consists of three main branches: implementing a QMC code and comparing it with existing software, and with results of previous studies, examining a minimal model of a TMD nanoribbon at the mean field level, applying our QMC code to simulate the model for TMD nanoribbons.

## Getting started

To run our Determinant QMC implementation simply open the following directory:


cd DETELECTRO-1.0/fast-version


and compile the code using


make


Then, you can run a simulation by running the command


./simulation 1 4 0 1 0 20384 512 128


which reproduces the results of the seminal paper "Discrete Hubbard-Stratonovich transformation for fermion lattice models", Phys Rev B, 28, 7, 1983, by J. E. Hirsch.

![alt-text][hirsch]

[hirsch]: https://github.com/fmonteir/qmc_master_code/blob/master/DETELECTRO-1.0/hirsch-reproduce.png

The arguments of *simulation* are described below.

To maximize the efficiency of the code, some simulation parameters must be known at compile time. Thus, to set them, you must provide arguments when running make. If you simply run make, you will obtain the default parameters that reproduce Hirsch's results for U = 4.


To change the number of sites, inverse Trotter error, inverse temperature, or the frequency of recomputing the Green's functions, type:


make clean


make nsites=\<Number of sites\> dt_inv=\<Inverse Trotter Error\> beta=\<Inverse Temperature\> green_afresh_freq=\<Frequency of Recomputing G\>


To run another simulation, simply type ./simulation followed by its arguments:


./simulation \<t\> \<U\> \<mu\> \<geom\> \<Ny\> \<Total Number of Sweeps (Space-Time)\> \<Number of Warm-up Sweeps (Space-Time)\>  \< Number of Auto-correlation Sweeps (Space-Time) \>

where the first argument is a hopping parameter, the second one is the on-site interaction, followed by the chemical potential. The next two parameters are related to the geometry of the model.
The last three parameters are related to the Monte Carlo method.


The program is prepared to handle a variety of geometries (listed below).
Input the number corresponding to the desired geometry:


(1)		  1D Periodic Chain

(2) 		1D Open Chain

(3) 		2D Periodic Square Lattice

(4) 		2D Open Square Lattice

(5) 		2D Periodic Rectangular Lattice

(6) 		2D Open Rectangular Lattice

(7) 		2D Periodic Triangular Lattice

(8) 		2D Nanoribbon Triangular Lattice

(9) 		2D Periodic Honeycomb Lattice

(10)		2D Honeycomb Nanoribbon

(11)		2D Honeycomb Hexagonal Dot

(12)		2D Honeycomb Triangular Dot

(13)		2D Honeycomb Rectangular Dot

(14)		2D Minimal model of a periodic TMD (MoS2) sample (Liu et al., Phys Rev B 88, 085433, 2013 ) - nsites includes orbital space, i.e. nsites=n_orbitals * n_spatial_sites.

(15)		2D Minimal model of a TMD (MoS2) nanoribbon (Liu et al., Phys Rev B 88, 085433, 2013 ) - nsites includes orbital space, i.e. nsites=n_orbitals * n_spatial_sites.

The Ny parameter is only meaningful for geometry options 5, 6, 8, 10, 13, and 15.

The results of the simulation will be saved in a directory named _temp-data_
and will be deleted once you run _make clean_. To save them enter the directory
_results_


cd results


and run


python save-simulation-data.py


In the _results/plot-src_ directory there are python scripts that use *matplotlib*
to plot the results of the simulations saved in _results/data_

For this work, we started by focusing on Molybdenum disulfide, MoS2.
This case of particular interest is obtained when one sets _geom=15_.
With this choice, we can explore how the properties of a MoS2 nanoribbon,
as obtained by a three-band minimal tight binding model are changed by adding
a Hubbard type interaction.

The hopping parameters allowing us to simulate a MoS2 nanoribbon are given in

[Liu2013](https://github.com/fmonteir/msc_references/blob/master/references/tmd/Liu2013.pdf)

When running _simulation_, we  should set t=-1, so as to normalize all the
hopping parameters to t0, i.e. setting t0=-1, and measuring everything in units of
t0.
Moreover, one must pay attention when setting the number of sites because it includes both
real and orbital space. For three orbitals, this means that if we set nsites=60,
we will be studying a 20-site system.

Other choices of the geom parameter allow the study of other TMDs (by changing the
  parameters of the minimal 3-band model introduced by Liu et al.)

### ALF-1.0, QUEST-1.4.9

Standard software used in the literature for benchmarking and comparison purposes.

### ANALYTICAL

Mean field studies, and results of basic analytical calculations for limiting cases of interest.

### DETELECTRO-1.0

Our own implementation of the determinant QMC algorithm to simulate the Hubbard model. The working name of the software is DETELECTRO. :)

In addition to our implementation of the BSS (Blakenbecler, Scalapino and Sugar) determinant QMC algorithm, this directory contains other interesting studies carried out in the context of this work, namely having to do with the sign problem, low temperature, and large size stabilization, etc.

## Built with

*C++*

*Python*

## Contributing

*Francisco Monteiro de Oliveira Brito* (CeFEMA, Instituto Superior Técnico de Lisboa, Centro de Física da Universidade do Porto)

## Versioning

v1 - latest update 14.08.2018

## Authors

*Francisco Monteiro de Oliveira Brito* (CeFEMA, Instituto Superior Técnico de Lisboa, Centro de Física da Universidade do Porto)

## Acknowledgements

*João M. V. P. Lopes, Eduardo V. Castro* - supervisors
