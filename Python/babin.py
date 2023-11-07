import sys
import re
import struct
import numpy as np

class bin :
    def __init__(self,fname, mode='rb') :
        """Create binary file object """
        self.n1=0
        self.d1=0
        self.o1=0
        self.n2=0
        self.d2=0
        self.o2=0
        self.file=open(fname,mode)
             
    def read(self,dim) :
        """Read binary data"""
        n=product(dim)
        data=np.fromfile(self.file, count=n ,dtype='float32')
        data=data.reshape(dim)
        self.file.close()
        return data 

    def readb(self) :
        """Read binary data"""
        data=np.fromfile(self.file, dtype='float32')
        self.file.close()
        return data 

    def write(self,data) :
        """Write binary data"""
        tmp = data.astype(np.float32)
        tmp.tofile(self.file)
        self.file.close()

def product(tuple1):
    """Calculates the product of a tuple"""
    prod = 1
    for x in tuple1:
        prod = prod * x
    return prod


