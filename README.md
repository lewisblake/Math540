# Parallel Scientific Computing

This repositiory contains sample code from a course titled Parallel Scientific Computing (Math 540) I took at Colorado School of Mines in Spring 2019.

The code is implemented in C++. There are four distinct projects in this repository.

## Assignment 1

This assignment executes serial code.

The files associated with this project are math540_hw1_1.cpp, math540_hw1_3a.cpp, math540_hw1_3b.cpp, functions.cpp, and functions.h.

## Assignment 2

The files associated with this project are assign2.cpp. This program is designed to utilize the C-implementation of the MPI library for parallel computations and the Eigen library for Linear Algebra.

The objective of this program is to evaluate a set of highly-oscialltory integrals. The PDF "Blake_Math540_HW2.pdf" acommpanies this program.

## Assignment 3

The files associated with this project are assign3.cpp This program is designed to utilize the C-implementation of the MPI library for parallel computations and the Eigen Library for Linear Algebra.

The objective of this program is to construct a large covariance matrix across processing cores. The PDF "lewisblake_assign3_math540.pdf" accompanies this program.

## Final Project

The files are located in the FinalProject folder, with the exeception of the Data folder, which is not included. In this folder there is an additional README.txt which outlines the functionality of the files. This code base builds the hierarchical domain decomposition of the multi-resolution approximation in C++. It was tested on the CSM HPC MIO.
