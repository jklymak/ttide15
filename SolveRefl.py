# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 15:50:08 2016

@author: jklymak
"""
import numpy as np
from vertmodes import  vertModes
import logging

def SolveRefl(k=0.,f =1.06e-4,omega=1.45e-4+1j*1e-6,wall=True,H=[],x=[],J=30,Nsq0=[],z0=[],
              Incoming=True,Forcing=[]):
    '''x,z,H,P=SolveRefl(k=0.)
     This one does variable Nsq
    ''' 
    log=logging.getLogger('SolveRefl')
    inv = np.linalg.inv
    ## Model Params:  
    k00=k
    dz = 1./(J+1)
    usebcs = True
    z = np.linspace(-1.,0.,J+1)
    dz = z[1]-z[0]
    #z = z[:-1]+0.5*dz
    zmid = z[:-1]+0.5*dz
    
    #lamsq = (np.real(omega)**2-f**2)/(Nsq0 - np.real(omega)**2)+0*z
    I = len(x)
    dx = x[2]-x[1]

    H[-2]=H[-1]; 
    H[1]=H[0]
    Hx = np.gradient(H,dx);Hxx = np.gradient(Hx,dx)
    # Just calculate Nsq once...
    
    Z0 = z[:,np.newaxis]*H[np.newaxis,:]
    Nsq = 0.*Z0
    for a in range(len(H)):
        Nsq[:,a]=np.interp(Z0[:,a],z0,Nsq0)
    dNsqdz=np.diff(Nsq,axis=0)/dz/H[np.newaxis,:]
    # now we want it on teh midpoints...
    Z0 = zmid[:,np.newaxis]*H[np.newaxis,:]
    Nsq = 0.*Z0
    for a in range(len(H)):
        Nsq[:,a]=np.interp(Z0[:,a],z0,Nsq0)
    
    N2 = np.interp(zmid*H[0],z0,Nsq0)
    psi,phi,ce,zpsi=vertModes(N2[::-1],dz*H[0])
    log.debug(ce[:10])
    #flip back because our matrices are from the bottom 
    psi = psi[::-1,:]*np.sqrt(dz*H[0])

    Z = np.diag(zmid,k=0)+0.*1j;eye = np.eye(J)+0.*1j
    ## To get started, we need P01, P02, and E and Ep...
    if 1:
        Amp=1.+1j*0.
        #pin1=pin1/np.max(np.abs(pin1))
        # Hmmm, OK this is hydrostatic.  What should it be?
        k1 = (np.real(omega)**2 - f**2)/ce[0]**2 -k00**2
        log.debug(np.pi*2./k1/1.e3,ce[:10])
        if k1<0:
            k1 = -1j*np.sqrt(-k1) # forcingd decays to left e^{kx}
        else:
            k1 = np.sqrt(k1) # wave is cominhg from right  e^{j kx}
        pin1 = Amp*psi[:,0]/psi[0,0]*np.exp(1j*k1*x[0])
        pin2 = pin1 + 1j*k1*dx
        log.debug(psi[0,0],pin1[0],pin2[0])
        pin2 = Amp*psi[:,0]/psi[0,0]*np.exp(1j*k1*x[1])  # = pin1 + dx *dPin/dx 
        log.debug(pin2[0])
    if not(Incoming):
        pin1=0.*pin1
        pin2=0.*pin2

    # LHS radiating boundnary condition.  Note that for k00 neq 0
    # k^2 can be negative, which implies a leftward decay...   
    nmodes = J-4
    E1 = psi[:,0:nmodes]+0.*1j
    K = np.zeros((nmodes,nmodes))*1j
    for j in range(nmodes):
        kk = (np.real(omega)**2 - f**2)/ce[j]**2 - k00**2
        if kk<0:
            kk = -1j*np.sqrt(-kk)  # decay to left e^(kx)
        else:
            kk = -np.sqrt(kk) # leftward! e^(-jkx)
        K[j,j] = 1j*kk*dx
    log.debug('lam',np.pi*2./K[0,0]/1e3*dx)
    alpha=[]
    beta=[]
    for k in range(I+1):
        alpha.append(np.zeros((J,J))*1j)    
        beta.append(np.zeros((J))*1j)
    log.debug(E1.dot(E1.transpose().conj()))
    ee = (E1.dot(K)).dot(E1.transpose().conj())
    alpha[0]=inv(eye+ee)
    beta[0]=pin1-(inv(eye+ee)).dot(pin2)

    P1 = np.zeros((J,J));
    P2 = np.zeros((J,J));
    for j in range(1,J-1):
        P1[j,j-1]=-1;P1[j,j+1]=1;
        P2[j,j-1]=1;P2[j,j]=-2;P2[j,j+1]=1;
    P1[0,1]=1;P1[J-1,J-2]=-1;
    P2[0,0]=-2; P2[0,1]=1;P2[J-1,J-2]=1;P2[J-1,J-1]=-2;
    ## So, that eliminates the error from the LHS....
    dxsq=dx**2
    for i in range(1,I-1):
        lamsq = (np.real(omega)**2-f**2)/(Nsq[:,i] - 0.*np.real(omega)**2)
        gamsq = lamsq/(Nsq[:,i] - 0.*np.real(omega)**2)
        lamsq = np.diag(lamsq,k=0)+1j*0.
        gamsq = np.diag(gamsq,k=0)+1j*0.
        G2 = -2*Hx[i]/H[i]*Z +0*1j      # G2  1/m
        G3 = -(Hxx[i]*H[i]-2*Hx[i]**2)/H[i]**2*Z + 0*1j  # 1/m^2
        G3 = G3 +  gamsq*dNsqdz[:,i]/H[i]  
        G4 = (Hx[i]**2*Z.dot(Z)-lamsq)/H[i]**2 +0*1j   # G4  1/m^2
        # get A, B , C , D
        
        D = np.zeros((J))*1j ## This needs to be set to something if you want internal forcing (versus an incoming wave)
        
        if not(Incoming):
            if len(Forcing)==0:
                if np.abs(x[i]-50e3)<1000000.e3:
                    D[2*J/4:3*J/4]+=1e-9/H[i]
            else:
                D = Forcing[:,i]
        A = eye*1./dxsq - G2.dot(P1)/4./dx/dz +0.*1j   #1/m^2
        B = -eye*(2./dxsq + k00**2) + G3.dot(P1)/2./dz +G4.dot(P2)/dz**2 +0.*1j   # 1/m^2 + 1/m^2 + 1/m^2
        C = eye*1./dxsq + G2.dot(P1)/4./dx/dz +0.*1j     #  1/m^2
        if True:          # use BC
            b1 = (lamsq[0,0] + zmid[0]*Hx[i]**2)/H[i]/H[i];  # 1/m^2 
            b2=-Hx[i]/H[i];  # 1/m 
            b3 = f*k00/omega*Hx[i]/H[i];  # 1/m^2

            # Note that this gives proper units to the BCs

            A[0,0]=-b2/2./dx; A[0,1]=0; # seafloor
            A[J-1,J-1]=0.;A[J-1,J-2]=0. # sea surface...
            
            B[0,0]= -b1/dz + b3;  B[0,1]= +b1/dz
            B[J-1,J-1]= lamsq[-1,-1]/H[i]**2/dz  
            B[J-1,J-2]=-lamsq[-1,-1]/H[i]**2/dz
             
            C[0,0]=b2/2./dx;C[0,1]=0;
            #C[0,0] += 1./dxsq;
            C[J-1,J-1]=0.;C[J-1,J-2]=0.;
            #C[J-1,J-1]+= 1./dxsq
        bb = inv(A.dot(alpha[i-1])+B)
        alpha[i]=-bb.dot(C)
        beta[i]=bb.dot(D)-(bb.dot(A.dot(beta[i-1])))

    # get P[I-1]
    P = np.zeros((J,I))+0.0*1j

    if wall:        
        log.debug("Wall")
        P[:,I-1] = inv(eye*(1.-f*k00*dx/omega)-alpha[I-2]).dot(beta[I-2])
    else: # radiating...
        N2 = np.interp(zmid*H[-1],z0,Nsq0)
        
        psi,phi,ce,zphi=vertModes(N2[::-1],dz*H[-1])
        psi = psi[::-1,:]*np.sqrt(dz*H[-1])

        E2 = psi[:,0:nmodes]+0.*1j
        K = np.zeros((nmodes,nmodes))*1j
        for j in range(nmodes):
            kk = (np.real(omega)**2 - f**2)/ce[j]**2 - k00**2
            if j==0:
                log.debug("kk[0]",kk)
            if kk<0:
                kk = 1j*np.sqrt(-kk)  # decay to right e^(-kx)
            else:
                kk = np.sqrt(kk) # rightward! e^(jkx)
            K[j,j] = 1j*kk*dx
            E2[:,j]=E2[:,j]*np.exp(1j*kk*x[-1])
        E2d=E2.dot(K)
        E2inv = E2.transpose().conj()
        P[:,I-1]= inv(eye-alpha[I-2]-E2d.dot(E2inv)).dot(beta[I-2])
    
    # back iterate to get P
    for i in range(I-2,-1,-1):
        P[:,i] = alpha[i].dot(P[:,i+1])+beta[i]

        
    return x,zmid,H,P,[E1,np.diag(K)]
   
