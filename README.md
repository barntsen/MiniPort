
# MiniPort 

MiniPort contains two simple mini applications written in C for
running on the Mi250x GPU. The mini applications solve the
acoustic wave equation using a simple finite-difference approach.
The purpose is to compare the performance of using managed versus
manual memory managment. 


The subdirectories contain the following:

## Hip-mmem

Source code for mini application using manual memory managment.
The source files contains:

 - ac2dmod.cpp : Main program
 - model.cpp   : Model definition
 - ac2d.cpp    : Solver routines
 - rec.cpp     : Recording data 
 - src.cpp     : Source definition
 
Script files for compiling and running:

 - mk.sh       : Compile the code
 - run.sh      : Create input data and execute
 - job.sh      : Submit job

## Hip-umem

Source code for mini application using unified memory managment.
The source files contains:

 - ac2dmod.cpp : Main program
 - model.cpp   : Model definition
 - ac2d.cpp    : Solver routines
 - rec.cpp     : Recording data 
 - src.cpp     : Source definition
 
Script files for compiling and running:

 - mk.sh       : Compile the code
 - run.sh      : Create input data and execute
 - job.sh      : Submit job

## Input and output data

The inputs to the mini applications are created by the run script
and consists of

- vp.bin  : Wave velocity model
- rho.vp  : Density model
- src.bin : Source pulse
- q.bin   : Q-values

The output are stored in two files: snp.bin and p.bin, where snp.bin
contains the complete acoustic wavefield, while p.bin only contain
the wavefield sampled at the receivers.

## Python
This directory contains python scripts for generating input data.
