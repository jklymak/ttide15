#!/usr/local/bin/python
# Filename: jmkdata.py

# from matplotlib import rc
from pylab import *
import numpy

def iwslope(om,f,Nsq):
    """
    def iwslope(om,f,Nsq)
    Return the slope of an internal wave ray (dz/dx) given f, Nsq and om.  
    """
    return  sqrt((om**2-f**2)/(Nsq-om**2))

