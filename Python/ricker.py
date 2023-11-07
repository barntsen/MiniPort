''' This script writes either a Gaussian source pulse, its first or
    second derivative (=ricker)  to the output file,
'''

import argparse
import numpy as np
import sys
from math import *
import babin as ba

#--------------------------------------------------------------------
#Get command line options
#---------------------------------------------------------------------
#Hack for negative floats as option
for i, arg in enumerate(sys.argv):
  if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg

# Read parameters from command line with parser
parser = argparse.ArgumentParser(description="Gaussian wavelets")
parser.add_argument("-nt",type=int, default=300,help="Number of time samples")
parser.add_argument("-dt",type=float, default=0.001,help="Time step")
parser.add_argument("-t0",type=float, default=0.1,help="Central time")
parser.add_argument("-f0",type=float, default=20.0,help="Central frequency")
parser.add_argument("-deriv",type=int, default=2,help="Derivative order")
parser.add_argument("-n",action='store_true',help='normalise')
parser.add_argument("fname",default="ricker.bin",help="Output binary file")
args  = parser.parse_args()


pi=3.14159
n     = args.n
nt    = args.nt
dt    = args.dt
w0    = 2.0*pi*args.f0
t0    = args.t0
deriv = args.deriv
fname = args.fname

pulse = np.zeros((nt))

for i in range(0,nt) :
    if(deriv == 0) :
        t = i*dt -t0
        arg = t*t*w0*w0/4.0
        pulse[i] = -exp(-arg)
    elif(deriv == 1): 
        t = i*dt -t0
        arg = t*t*w0*w0/4.0
        arg2 = 2*t*w0*w0/4.0
        pulse[i] = arg2*exp(-arg)
    elif(deriv == 2): 
        t = i*dt -t0
        arg = t*t*w0*w0/4.0
        arg2 = t*t*w0*w0/2.0
        if(n is False):
            pulse[i] = (1-arg2)*exp(-arg)
        else:
            pulse[i] = 2.0*w0*w0*(1-arg2)*exp(-arg)
 

fp = ba.bin(fname,"wb")
fp.write(pulse)
