#!/bin/sh

# Test script for PyAc2d. 

export HSA_XNACK=1
export HSA_ENABLE_SDMA=0

BIN=../Python

#Create wavelet
nt=1501 #No of samples
$BIN/ricker -nt $nt -f0 30.0 -t0 0.100 -dt 0.0005 src.bin 

n1=256
n2=256
#Create vp
$BIN/spike -n1 $n1 -n2 $n2 -val 2500.0 vp.bin

#Create rho 
$BIN/spike -n1 $n1 -n2 $n2 -val 1000.0 rho.bin

#Create Q 
$BIN/spike -n1 $n1 -n2 $n2 -val 100000.0 q.bin

#Run modelling
./ac2dmod > log.txt


