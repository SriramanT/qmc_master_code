#  What this code is meant to do

Here I try to reproduce the implementation discussed by Bai et al.

The code in this folder represents the completely naive implementation. I implemented an efficient update of the Green's functions that shall serve to perform measurements. A prominent problem is that this update involves wrapping G-matrices with B-matrices and this is numerically unstable if done in the naive way as in this code.

The plots clearly show that the convergence is happening in both the naive, brute force case and in the more efficient approach. Moreover, we can clearly see that numerical instability arising after just a few MC sweeps (NxL HS-field).

# What is missing?


This code lacks the implementation of the measurements (which is not very hard to make, since they depend exclusively on the Green's function). It also lacks a calculation of warm-up and auto-correlation times.

Something that is not yet clear to me, though, is whether the G-matrices or the G-hat-matrices are used when measuring within the efficient approach.

## Reproducibility

Care must taken when making plots so that the results are reproducible. This should be built into the code and the directory structure of the folder.

## Object-oriented

It would be nice to define classes and to give the code a better structure.

# Why does v1 stop here?

I felt like since I arrived at a point were things are starting to come together and make sense, I should document the progress by taking a step back to reflect upon what should be improved, changed and where to go from here.

