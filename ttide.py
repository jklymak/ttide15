
# coding: utf-8

# In[1]:

import numpy as np
import matfile as mf
import hdf5storage
import cPickle as pickle
from multiprocessing import Pool
import itertools,time

# get the exact shelves:



def SolveRefl(k=0.,f =1.06e-4,omega=1.45e-4+1j*1e-6,Nsq0=0.01**2,wall=True,H=[],x=[],J=30):
    '''x,z,H,P=SolveRefl(k=0.)
     
    ''' 
    inv = np.linalg.inv
    ## Model Params:  
    k00=k
    dz = 1./(J+1)
    usebcs = True
    if usebcs:
        shape=np.cos
        z = np.linspace(-1.,0.,J+1)
        z = z[:-1]+(z[1]-z[0])/2.
        
    else:
        shape = np.sin
        z = np.arange(-dz*(J),-dz+0.00001,dz)
        
    dz = z[1]-z[0]
    lamsq = (omega**2-f**2)/(Nsq0 - omega**2)+0*z

    if len(x)==0:
        I=64
        x = np.linspace(0.,60.1e3,I)
    else:
        I = len(x)
    dx = x[2]-x[1]
    # Some topo:
    if len(H)==0:
        H=3200-(3200-120)/2.*np.tanh(0.12e-3*(x-30e3))-(3200-120)/2;
    #H=225.0 - 50.0*np.exp(-((x-40.e3)/3.e3)**2)
    #H = 225+0*H
    H[-2]=H[-1]; 
    H[1]=H[0]

    Hx = np.gradient(H,dx);Hxx = np.gradient(Hx,dx)
    
    # set up matrix multipliers...
    P1 = np.zeros((J,J));
    P2 = np.zeros((J,J));
    for j in range(1,J-1):
        P1[j,j-1]=-1;P1[j,j+1]=1;
        P2[j,j-1]=1;P2[j,j]=-2;P2[j,j+1]=1;
    P1[0,1]=1;P1[J-1,J-2]=-1;
    P2[0,0]=-2; P2[0,1]=1;P2[J-1,J-2]=1;P2[J-1,J-1]=-2;

    Z = np.diag(z,k=0)+0.*1j;eye = np.eye(J)+0.*1j
    ## To get started, we need P01, P02, and E and Ep...

    ### P0:  Amp*exp(j(kn x - om t + k y))phi(z)
    Amp = 0.+0.*1j
    pin1 = Amp*shape(z*np.pi)
    k1 = (np.pi/H[0])**2 * lamsq[0] - k00**2
    if k1<0:
        k1 = -1j*np.sqrt(-k1) # forcingd decays to left
    else:
        k1 = np.sqrt(k1) # wave is cominhg from right
    pin2 = pin1*(1+1j*k1*dx)  # = pin1 + dx *dPin/dx 
    
    ### E and Ep
    # for general Nsq we need to solve the eignevalue problem.  Lets be lazy here.
    nmodes = J-2
    E1 = np.zeros((J,nmodes))*1j
    K = np.zeros((nmodes,nmodes))*1j
    for j in range(nmodes):
        E1[:,j]=-np.sqrt(2.*dz)*shape((j+1)*z*np.pi)
        kk = ((j+1)*np.pi/H[0])**2 * lamsq[0] - k00**2 + 0*1j
        kk = np.sqrt(kk)
        ## Note that kk can be complex:
        if np.real(kk)>0:
            # choose the negative so that we are radiating offshore:
            kk = -kk
        if 0:
            if kk<0:
                kk=-1j*np.sqrt(-kk)  # decay to left
            else:
                kk = -np.sqrt(kk) # leftward!
        K[j,j] = 1j*kk*dx
    alpha=[]
    beta=[]
    for k in range(I+1):
        alpha.append(np.zeros((J,J))*1j)    
        beta.append(np.zeros((J))*1j)
    ee = E1.dot(K).dot(E1.transpose().conj())
    alpha[0]=inv(eye+ee)
    beta[0]=pin1-inv(eye+ee).dot(pin2)

    ## So, that eliminates the error from the LHS....
    gamsq = lamsq/(Nsq0 - omega**2)
    dxsq=dx**2
    for i in range(1,I-1):
        #print i
        G2 = -2*Hx[i]/H[i]*Z                               # G2 
        G3 = -(Hxx[i]*H[i]-2*Hx[i]**2)/H[i]**2*Z
        G4 = (Hx[i]**2*Z.dot(Z)-np.diag(lamsq))/H[i]**2   # G4
        # get A, B , C , D
        D = np.zeros((J))*1j## This needs to be set to something if you want internal forcing (versus an incoming wave)
        if np.abs(x[i]-10e3)<0.5e3:
            D[:J/3]+=1./1000.
            D[J/3:]+=0./1000.
                
        A = eye*1./dxsq - G2.dot(P1)/4./dx/dz
        B = -eye*(2./dxsq + k00**2) + G3.dot(P1)/2./dz +G4.dot(P2)/dz**2
        C = eye*1./dxsq + G2.dot(P1)/4./dx/dz
        if usebcs:
            b1 = (lamsq[0]-Hx[i]**2)/H[i]; b2=Hx[i]; b3 = -f*k00/omega*Hx[i];
            
            A[0,0]=b2/2./dx; A[0,1]=0; # seafloor
            A[J-1,J-1]=1./4./dx**2;A[J-1,J-2]=0. # sea surface...
            #A[0,0]=1./4./dx**2;A[0,1]=0. # sea floor...
            
            B[0,0]= -b1/dz - b3;B[0,1]=+b1/dz
            B[J-1,J-1]=(-2./4./dx**2-k00**2) + lamsq[0]/dz;  B[J-1,J-2]=-lamsq[0]/dz
            #B[0,0]=(-2./4./dx**2-k00**2) + lamsq[0]/dz;  B[0,1]=-lamsq[0]/dz
            
            C[0,0]=-b2/2./dx;C[0,1]=0;
            C[J-1,J-1]=1./4./dx**2;C[J-1,J-2]=0.;
            #C[0,0]=1./4./dx**2;C[0,1]=0.;
        
        bb = np.linalg.inv(A.dot(alpha[i-1])+B)
        alpha[i]=-bb.dot(C)
        beta[i]=bb.dot(D)-(bb.dot(A).dot(beta[i-1]))

    # get P[I-1]
    P = np.zeros((J,I))+0.0*1j
    inv = np.linalg.inv

    if wall:        
        P[:,I-1] = inv(eye*(1.-f*k00*dx/omega)-alpha[I-2]).dot(beta[I-2])
    else: # radiating...
        E2 = np.zeros((J,nmodes))*1j
        E2d = np.zeros((J,nmodes))*1j
        K = np.zeros((nmodes,nmodes))*1j
        for j in range(nmodes):
            kk = ((j+1)*np.pi/H[I-1])**2 * lamsq[0] - k00**2
            kk = np.sqrt(kk)
            
            E2[:,j]=-np.sqrt(2.*dz)*shape((j+1)*z*np.pi)*np.exp(1j*kk*x[-1])
            K[j,j]= 1j*kk*dx
            E2d[:,j]=1j*kk*E2[:,j]
        E2d=E2.dot(K)
        E2inv = E2.transpose().conj()
        #phi(:,I)=inv(eye(J)-Hx(I)/H(I)/2/dz*dx*Z*P1-alp(:,:,I-1)-dx*E2d*inv(E2))*bta(:,I-1);
        P[:,I-1]= inv(eye-alpha[I-2]-E2d.dot(E2inv)).dot(beta[I-2])
        #P[0,I-1]=P[1,I-1]
        #P[-1,I-1]=P[-2,I-1]
        
    for i in range(I-2,-1,-1):
        P[:,i] = alpha[i].dot(P[:,i+1])+beta[i]
    return x,z,H,P


# In[12]:



def getRes(Hk):
    
    H=Hk[0]
    k=Hk[1]
    N0=Hk[2]
    xmod=Hk[3]
    f=1.e-4
    omega = np.pi*2./12.4/3600.

    x,z,H,P=SolveRefl(k=k,Nsq0=N0**2,omega=omega+1j*1e-6,f=f,wall=True,x=xmod*1e3,H=H,J=32*2)
    res=1.e-99
    for num,dt in enumerate(np.linspace(0.,np.pi,100)):
        ind = np.where(x>90e3)[0]
        re=np.mean(np.abs(np.real(P[:,ind]*np.exp(1j*dt)))*H[np.newaxis,ind])
        if re>res:
            res=re
    return res,k

if __name__ == '__main__':
    todo = ['Shelf100km','Shelf1km03','Shelf1km04','Shelf020km']
    Hwkb=[]

    xmod=np.linspace(0.,170.,140)
    for nn,td in enumerate(todo):
        D = hdf5storage.loadmat('../ttide15/Tas3d/%s/Diags0360.mat'%td)
        H = D['Depth'].transpose()[:,-1]
        x = -D['x']/1e3-30.+170
        ray=mf.loadmatbunch ('../TasmaniaRays.mat')
        ray = ray['ray']
        z=np.arange(0.,5000.,2.)
        
        N=  np.interp(z,ray['z'],np.sqrt(ray['N2']))
        N0 = np.mean(N)
        zwkb=np.cumsum(N/N0*2)
        hh = np.interp(H,z,zwkb)
        Hwkb.append(np.interp(-xmod,-x,hh))
        
    lam = np.linspace(30,350.0,100)*1e3
    ks = np.pi*2./lam
    kims = np.linspace(0,np.pi*2./40e3,2)
    
    resp=[]
    ## Scan im/re k space:
    hin=2

    p = Pool(1)
    for hin in range(1):
        resp.append(np.zeros((len(kims),len(ks))))
        t0 = time.time()
        print hin
        for m,kim in enumerate(kims):
            print m
            ## try pool
            k = ks+1j*kim # propagate N and decay to north...
            Hk= itertools.izip(itertools.repeat(Hwkb[hin]),k,itertools.repeat(N0),itertools.repeat(xmod))
            r = p.map(getRes,Hk)
            rr = np.zeros(len(r))
            kk = np.zeros(len(r))*1j
            for nnn in range(len(r)):
                rr[nnn]=r[nnn][0]
                kk[nnn]=r[nnn][1]


            
            resp[hin][m,:]=rr[np.argsort(np.real(kk))]

        t1=time.time()
        print(t1-t0)

    pickle.dump(resp,open('resp.pickle','wb'))

