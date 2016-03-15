
# coding: utf-8

# In[1]:

import numpy as np
import matfile as mf
import hdf5storage
import cPickle as pickle
from multiprocessing import Pool
import itertools,time
import sys

# get the exact shelves:

def getRes(Hk):

    import SolveRefl as sr

    H=Hk[0]
    k=Hk[1]
    xmod=Hk[2]
    f=1.e-4
    omega = Hk[3]

    x,z,H,P,debug=sr.SolveRefl(k=k,Nsq0=Nsq0,z0=z0,omega=omega*f,f=f,wall=True,x=xmod*1e3,
                           H=H,J=192,Incoming=False)
    ind = np.where(x>(x[-1]-70e3))[0]
    res=np.mean(np.abs(np.real(P[:,ind]*np.exp(1j*0.)))*H[np.newaxis,ind])
    return res,k

if __name__ == '__main__':
    todo = ['Shelf100km','Shelf1km03','Shelf1km04','Shelf020km']
    Hwkb=[]

    xmod=np.linspace(0.,260.,261)
    for nn,td in enumerate(todo):
        D = hdf5storage.loadmat('../ttide15/Tas3d/%s/Diags0360.mat'%td)
        H = D['Depth'].transpose()[:,-1]
        x = -D['x']/1e3-30.+xmod[-1]
        ray=mf.loadmatbunch ('../TasmaniaRays.mat')
        ray = ray['ray']
        z=np.arange(0.,5000.,2.)

        N=  np.interp(z,ray['z'],np.sqrt(ray['N2']))
        N0 = np.mean(N)
        Hwkb.append(np.interp(-xmod,-x,H))

    lam = np.linspace(30,180.0,140)*1e3
    ks = np.pi*2./lam
    oms = np.linspace(1.,2.,3)
    oms = oms[1:]
    z0 = np.linspace(-6500,0.,1000)
    # Nsq0 = np.interp(z0,-(ray['z'][2:][::-1]*1.),ray['N2'][2:][::-1])
    #Nsq0 = np.interp(z0,-(ray['z'][2:][::-1]*1.),np.sort(ray['N2'][2:][::-1]))
    Nsq0 = 2e-5*np.exp(z0/1000.)

    resp=[]
    ## Scan im/re k space:
    hin=2

    p = Pool(3)
    for hin in range(4):
        
        resp.append(np.zeros((len(oms),len(ks))))
        t0 = time.time()
        print hin
        for m,om in enumerate(oms):
            sys.stdout.write('%d'%m)
            ## try pool
            k = ks # propagate N and decay to north...
            
            Hk= itertools.izip(itertools.repeat(Hwkb[hin]),k,itertools.repeat(xmod),itertools.repeat(om))
            r = p.map(getRes,Hk)
            rr = np.zeros(len(r))
            kk = np.zeros(len(r))*1j
            for nnn in range(len(r)):
                rr[nnn]=r[nnn][0]
                kk[nnn]=r[nnn][1]


            
            resp[hin][m,:]=rr[np.argsort(np.real(kk))]
            sys.stdout.write('\r')

        t1=time.time()
        print(t1-t0)

        pickle.dump(resp,open('respkomTest.pickle','wb'))

