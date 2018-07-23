#  What this code is meant to do

Here I try to reproduce the implementation discussed by Bai et al. in the following document:

https://github.com/fmonteir/qmc_master/blob/master/reading/finite-temperature-afqmc/Bai2006.pdf

The code in this folder corresponds to the completely naive implementation. An efficient update of the Green's functions that shall serve to perform measurements was also implemented. A prominent problem is that this update involves wrapping G-matrices with B-matrices and this is numerically unstable if done in the naive way as in this code (although the cost is the same: O(N^3) ).

The plots clearly show that the convergence is happening in both the naive, brute force case and in the more efficient approach. Moreover, we can clearly see that numerical instability arising after **just a few** MC sweeps (NxL HS-field).

# Compiling the code

Using g++, on my machine (a 2015 MacBook Pro), it suffices to run

*g++ -std=c++11 -o hubbard-1d-finiteT -I includes hubbard-1d-finiteT.cpp generateHSmatrix.cpp createHoppingMatrix.cpp prints.cpp matrices.cpp*

I placed the linear algebra library that is used in the code in the default folder where C++ looks for includes when compiling any code. That is why there is no explicit reference to it in the command above.

# Making plots

The *plots* folder has 3 .txt files corresponding to the data of the previously run simulation. There are two python scripts (in both .py and .ipynb formats). One makes the plot corresponding to the previous simulation. The other reproduces a previously saved plot. There are already two example folders with plots and the corresponding data.

# Outputs

The code outputs a list of the parameters used in the simulation, and two lists of statistical weights of the accepted HS configurations: one is obtained by naively computing the B and M-matrices throughout the simulation, the other is obtained by performing the more efficient updates. The output is saved in then saved in the folder *plots*.

Running the script plots.ipynb, we create a folder with a name that contains the parameters of the simulation. In this folder, the data that is used to make the plots is saved (the weights), and a plot is created. Depending on the number of MC sweeps, we either see agreement between the blue curve (representing the naive weights), and the blue one (representing the ones computed via update), or we the naive weights in blue converging and the ones obtained by the updates running off.


# What is missing?

This code lacks the implementation of the measurements (which is not very hard to make, since they depend exclusively on the Green's function). It also lacks a calculation of warm-up and auto-correlation times.

Something that is not yet clear to me, though, is whether the G-matrices or the G-hat-matrices are used when measuring within the efficient approach.

## Reproducibility

Care must taken when making plots so that the results are reproducible. This should be built into the code and the directory structure of the folder.

## Object-oriented

It would be nice to define classes and to give the code a better structure.

# Why does v1 stop here?

I felt like since I arrived at a point were things are starting to come together and make sense, I should document the progress by taking a step back to reflect upon what should be improved, changed and where to go from here.



