#!/bin/sh

# Test script for PyAc2d. 

export HSA_XNACK=1
export HSA_ENABLE_SDMA=0


#Create wavelet
nt=1501 #No of samples
ricker -nt $nt -f0 30.0 -t0 0.100 -dt 0.0005 src.bin 

n1=256
n2=256
#Create vp
spike -n1 $n1 -n2 $n2 -val 2500.0 vp.bin

#Create rho 
spike -n1 $n1 -n2 $n2 -val 1000.0 rho.bin

#Create Q 
spike -n1 $n1 -n2 $n2 -val 100000.0 q.bin

#Run modelling
./ac2dmod > log.txt

#../../Scripts/snp.sh

