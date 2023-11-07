""" Spike : generate spikes in a 3D binary cube """

from __future__ import print_function
import sys
import argparse
import numpy as np
import babin as ba

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# Define options

#---------------------------------------------------------------------
#Hack for negative floats as option
for i, arg in enumerate(sys.argv):
  if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg
#---------------------------------------------------------------------
parser = argparse.ArgumentParser(description="program to window a binary cube")
parser.add_argument("fout",help="Output binary file")
parser.add_argument("-n1",dest="n1",type=int,help="first dimension of data")
parser.add_argument("-n2",dest="n2",type=int,default=1,help="second dimension of data")
parser.add_argument("-n3",dest="n3",type=int,default=1,help="third dimension of data")
parser.add_argument("-f1",dest="f1",type=int,default=0,help="start of spike in first dimension")
parser.add_argument("-f2",dest="f2",type=int,default=0,help="start of spike in second dimension")
parser.add_argument("-f3",dest="f3",type=int,default=0,help="start of spike in third dimension")
parser.add_argument("-l1",dest="l1",type=int,help="length of spike in first dimension")
parser.add_argument("-l2",dest="l2",type=int,help="length of spike in second dimension")
parser.add_argument("-l3",dest="l3",type=int,help="length of spike in third dimension")
parser.add_argument("-i1",dest="i1",type=int,help="increment in first dimension")
parser.add_argument("-i2",dest="i2",type=int,help="increment in second dimension")
parser.add_argument("-i3",dest="i3",type=int,help="increment in third dimension")
parser.add_argument("-val",dest="val",type=float,default=1.0,help="value of spike")

#Parse arguments
args = parser.parse_args()

# Axis
if args.n1 is not None:
    n1 = args.n1
else :
    sys.exit("Missing n1!")

if args.l1 is not None:
    l1 = args.l1
else :
    l1 = args.n1

if args.i1 is not None:
    i1 = args.i1
else :
    i1 = 1

if args.l2 is not None:
    l2 = args.l2
else :
    l2 = args.n2

if args.i2 is not None:
    i2 = args.i2
else :
    i2 = 1

if args.l3 is not None:
    l3 = args.l3
else :
    l3 = args.n3

if args.i3 is not None:
    i3 = args.i2
else :
    i3 = 1

n2 = args.n2
n3 = args.n3
f1 = args.f1
f2 = args.f2
f3 = args.f3

if n1 <= 0 :
    sys.exit("n1 0 or negative \n")
if n2 <= 0 :
    sys.exit("n2 0 or negative \n")
if n3 <= 0 :
    sys.exit("n3 0 or negative \n")

if(f1 < 0):
    sys.exit("f1  negative \n")
if(f2 < 0):
    sys.exit("f2  negative \n")
if(f3 < 0):
    sys.exit("f3  negative \n")

if(l1 <= 0):
    sys.exit("l1 0 or negative\n")
if(l2 <= 0):
    sys.exit("l2 0 or negative\n")
if(l3 <= 0):
    sys.exit("l3 0 or negative\n")

if(f1+l1 > n1):
    sys.exit("f1+l1 too large\n")
if(f2+l2 > n2):
    sys.exit("f2+l2 too large\n")
if(f3+l3 > n3):
    sys.exit("f3+l3 too large\n")

#Create the data
data=np.zeros((n3,n2,n1))
print(data.shape)

#Insert spike
data[f3:f3+l3:i3,f2:f2+l2:i2,f1:f1+l1:i3] = args.val

#write the data
fout = ba.bin(args.fout,"w")
fout.write(data)

