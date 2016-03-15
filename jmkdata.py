#!/usr/local/bin/python
# Filename: jmkdata.py

# from matplotlib import rc
from pylab import *
import numpy
import numpy as np

def gappyinterp(tnew,t,x,maxdt):
    """
    gappyinterp(tnew,t,x,maxdt):
    if a gap is greater than dt, NaNs are returned in the data gap.
    """
    dt=diff(t)
    bad=squeeze(where(abs(dt)>maxdt))
    xnew = ma.array(interp(tnew,t,x,left=NaN,right=NaN))
    ind = where((tnew<t[0]) | (tnew>t[-1]))
    xnew[ind]=ma.masked
    if bad.size>1:
        for nn in range(0,shape(bad)[0],1):
            b=bad[nn]
            inds = where((tnew>t[b]) & (tnew<t[b+1]))
            xnew[inds]=ma.masked
    return xnew 

def bindata(binx,biny,x,y,X):
    """
    bindata(binx,biny,x,y,X):
    returns data Xnew binned into binx and biny bins
    meanX,varX,nX=bindata(binx,biny,x,y,X)
    """

    import numpy as np

    M = size(biny,0)-1
    N = size(binx,0)-1

    good = ~np.isnan(x+y+X)
    x=x[good]
    y=y[good]
    X=X[good]
    
    meanX = np.zeros((M,N))
    varX = np.zeros((M,N))
    nX = np.zeros((M,N))
    
    #print shape(binx)
    #print shape(arange(0,N,1.0))
    
    indx = floor(np.interp(x,binx,arange(0,N+1,1.0)))
    indx[indx<0]=0
    indx[indx>N-1]=N-1

    indy = floor(np.interp(y,biny,arange(0,M+1,1.0)))
    #print shape(X)
    indy[indy<0]=0
    indy[indy>M-1]=M-1
    for ind in range(len(X)):
        meanX[indy[ind],indx[ind]]+=X[ind]
        nX[indy[ind],indx[ind]]+=1
    nX[nX==0.0]=NaN;
    meanX = meanX/nX;

    for ind in range(len(x)):
        varX[indy[ind],indx[ind]]+=(meanX[indy[ind],indx[ind]]-X[ind])**2
    varX=varX/nX
    # mask
    meanX=np.ma.masked_where(isnan(nX),meanX)
    nX=np.ma.masked_where(isnan(nX),nX)
    varX=np.ma.masked_where(isnan(nX),varX)
    
    
    return meanX,varX,nX

def gappy_fill(x):
    """ gappy_fill(x) fills gaps in x"""
    xnew=x.copy()
    for ind in range(shape(x)[0]):
        good = where(~isnan(x[ind,:]))[0]
        if len(good)>1:    
            toin=range(good[0],good[-1])
            xnew[ind,toin]=interp(toin,good,x[ind,good])
    return xnew

def vertModesOld(N2,dz): 
    """" psi,zpsi,phi,ce=vertModesOld(N2,dz)
    
    vertModes moved to vertmodes.VertModes

    Compute the vertical eigen modes of the internal wave solution on a flat bottom
    
    Inputs: N2 is buoyancy frequency squared (rad^2/s^2) as an 1-D array. 
    
            dz is a single value, and the distance (in meters) between
            the N2 estimates
            
            If there are M values of N2, the first one is assumed to
            be at dz/2 deep, and the last one is H-dz/2 deep.  The
            water column is assumed to be H=M*dz deep.  No gaps are
            allowed, and N2>0 everywhere.
            
            Note that for M>200 or so, this gets quite slow, and you
            may ask why you are fitting more than 200 vertical modes.
            
    Outputs: psi (M,M) is the horizontal structure function at
             z=dz/2,3dz/2...H-dz/2.  psi is normalized so that
             sum(psi^2 dz) = 1.  If you interpolate psi onto a
             different value of dz, you probably want to renormalize.
             For internal waves, psi is approriate for velocity and
             pressure vertical structure.
             
             zpsi (M) is the locations of the psi values in the
             vertical
             
             phi (M+1,M) is the vertical integral of psi (phi = int psi dz) and represents 
             the vertical velocity structure at z = 0,dz,2dz...H
             
             ce is the non-rotating phase speed of the waves in m/s.
             
    Notes: This solves 
        d/dz(1/N**2 d psi/dz) + (1/ce**2)psi_{zz} = 0
    subject to a boundary condition of no derivative at the boundaries.  
    psi_{z}(0)=0 (rigid lid approx)
    psi_{z}(H)=0
    It is solved as an eigenvalue problem.  
             
    J. Klymak (Based on code by Gabriel Vecchi)           
    """
    if size(dz)>1:
        error('dz must be a constant')
    M = shape(N2)[0]
    # get psi:
    D = zeros((M,M))
    # surface:
    surffac = 1.
    print "surffac=%f"%surffac
    D[0,0] = -surffac/N2[0]
    D[0,1] = surffac/N2[0]
    # interior:
    for i in arange(1,M-1): 
        D[i,i-1] = 1./N2[i-1]
        D[i,i] = -1./N2[i-1]-1./N2[i]
        D[i,i+1] = 1./N2[i]
    # bottom:
    D[-1,-2] = surffac/N2[-1]
    D[-1,-1] = -surffac/N2[-1]
    D=-D/dz/dz
    ce,psi = numpy.linalg.eig(D)
    # psi is such that sum(psi^2)=1 but we want sum(psi^2 dz)=1.
    psi = psi/sqrt(dz)
    ce = 1./sqrt(ce)
    ind=argsort(-ce)
    if ce[ind][0]>50.:
        ind = ind[1:]
    ce = ce[ind]
    psi = psi[:,ind]
    zpsi = np.arange(dz/2.,M*dz,dz)
    #print shape(zpsi)
    INT = tril(ones((M,M)));
    INT = INT - 0.5*(eye(M));
    INT[:,0] = INT[:,0] - 0.5;
    INT = INT*dz
    phi = np.zeros((M+1,len(ind)))
    phi[1:-1,:] = np.diff(psi,axis=0)
    phi = phi/np.sqrt(np.sum(phi**2,axis=0)*dz)
    zphi = np.arange(0.,dz*(M+1),dz)
    #phi=phi[:,ind[0:-2]]
    
    return psi,zpsi,phi,zphi,ce
