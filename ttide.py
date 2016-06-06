# plotting the energy budget from a structure D.
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import scipy.signal as signal
def plotEnergyBudget(D):
    xl = [-40.,100.]
    yl=[0.,400.]

    Pu=D['uPbc']
    Pv = D['vPbc']
    H=D['Depth']
    x=D['x']/1e3
    y=D['y']/1e3
    
    cmap=cm.get_cmap('RdBu_r')

    divPbc = np.diff(Pu[:-1,:],axis=1)/np.diff(x)
    divPbc+=(np.diff(Pv[:,:-1],axis=0).T/np.diff(y)).T
    #divPbc=divPbc/1000.

    dEdt = 1000*(D['Ebc']-D['Ebc0'])/12.4/3600.

    #    djmkfigure(2,0.5)
    fig = plt.figure(figsize=(10,5.3))
    gs=gridspec.GridSpec(1,4,right=0.87,left=0.085,bottom=0.08,wspace=0.07,top=0.97)
      
    kernel=np.ones((4,4))/16.
    ax=[0,0,0,0]
    Z=[0,0,0,0]
    Z[0]=signal.convolve2d(D['Conv']*1000.*1000.,kernel,mode='same')
    Z[1]=-1000*signal.convolve2d(divPbc,kernel,mode='same')
    Z[2]=-1000*signal.convolve2d(dEdt,kernel,mode='same')
    Z[3]=-1000.*signal.convolve2d(divPbc-1000.*D['Conv'][:-1,:-1]+dEdt[:-1,:-1],kernel,mode='same')
    tit=['BT-BC conv.','BC Convergence','dE/dt','Diss. $[mW/m^2]$']

    for nn in range(4):
        ax[nn]=plt.subplot(gs[nn])
    
        pcm=ax[nn].pcolormesh(x,y,Z[nn],
                     rasterized=True,cmap=cmap)
        pcm.set_clim(np.array([-1.,1.])/20./4.*1000.)
        plotdepthcont(ax[nn],x,y,H)
        plotBox(ax[nn],x,y,xl,yl)
        ax[nn].set_aspect(1.)
        plt.title(tit[nn],fontsize=12)
        plt.xlabel('X [km]')
        plt.ylabel('Y [km]')
   
    colorbarRight(pcm,gs,fig,width=0.012,shrink=0.55,extend='both')
    for aa in ax:
        print aa
        aa.set_ylim([-80.,500.])
        aa.set_xticks(np.arange(-200.,200.,100.))
        aa.set_xlim([-99.,180])
        
    for aa in ax[1:]:
        aa.set_yticklabels('')
        aa.set_ylabel('')
    sublabel(np.array(ax),fontsize=12)
    return fig

def plotBox(ax,x,y,xl,yl):
    inx=np.where((x>xl[0]) & (x<=xl[1]))[0]
    iny=np.where((y>yl[0]) & (y<=yl[1]))[0]
    ax.plot(x[inx[[0,-1,-1,0,0]]],y[iny[[0,0,-1,-1,0]]],color='g',linewidth=1.)
def plotsponge():
    spongew=40
    plot(x[spongew-1]*array([1,1]),y[[0,-1]],'c',linewidth=1,alpha=0.5)
    plot(x[-spongew]*array([1,1]),y[[0,-1]],'c',linewidth=1,alpha=0.5)
    plot(x[[0,-1]],y[spongew-1]*array([1,1]),'c',linewidth=1,alpha=0.5)
    plot(x[[0,-1]],y[-spongew]*array([1,1]),'c',linewidth=1,alpha=0.5)
def plotdepthcont(ax,x,y,H):
    ax.contour(x,y,-H,[-250.,-3000,-2000,-1000,-4000,-10000],colors='k',linestyles='solid',alpha=0.5,linewidth=1.5)
    ax.contourf(x,y,-H,[-1.,0.],colors=[[0.2,0.45,0.2]],linestyles='solid',alpha=1.,linewidth=1.5)
    #inx = where(diff(x)<=1.)[0]
    #iny = where(diff(y)<=1.)[0]
    #plot(x[inx[[0,-1,-1,0,0]]],y[iny[[0,0,-1,-1,0]]],color='g',linewidth=1.)
#pcm=pcolormeshRdBu(x,y,P)

def sublabel(axs,fontsize=9):
    '''
    sublabel(axs,fontsize=9):
    '''
    for nn,ax in enumerate(axs.flatten()):
        ax.text(0.05,1.-0.07,'%c)'%chr(ord('a')+nn),
                fontsize=fontsize,transform = ax.transAxes,
                color='#555555',
               bbox=dict(facecolor='w', edgecolor='None',
                        alpha=0.85))


# LatLon to Model
def lonlat2modxy(lon,lat):
    
    Lat0=-44.
    Lon0=148
    kmpernm = 1.8532

    x=(lon-Lon0)*kmpernm*60.*np.cos(Lat0*np.pi/180)
    y=(lat-Lat0)*kmpernm*60.
    xx=x+1j*y
    xx=xx*np.exp(1j*12.*np.pi/180.)
    return np.real(xx),np.imag(xx)
def modxy2lonlat(x,y):
    Lat0=-44.
    Lon0=148
    kmpernm = 1.8532
    xx=x+1j*y
    xx=xx*np.exp(-1j*12.*np.pi/180.)
    lon = np.real(xx)/kmpernm/60/np.cos(Lat0*np.pi/180.)+Lon0
    lat = np.imag(xx)/kmpernm/60+Lat0
    return lon,lat

def colorbarRight(pcm,ax,fig,shrink=0.7,width=0.025,gap=0.03,**kwargs):
    '''
    def colorbarRight(pcm,ax,fig,shrink=0.7,width=0.05,gap=0.02)
    
    Position colorbar to the right of axis 'ax' with colors from artist pcm.
    ax can be an array of axes such as that returned by "subplots".
    
    ax can also be a GridSpec, in which case the colorbar is centered to the
    right of the grid.  
    
    Defaults might no leave enough room for the colorbar on the right side, so 
    you should probably use subplots_adjust() or gridspec_update() to make more 
    space to the right:
    
    # with subplots:
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(2,2)
    fig.subplots_adjust(right=0.87)
    for axx in ax.flatten():
        pcm=axx.pcolormesh(rand(10,10))
    colorbarRight(pcm,ax,fig,extend='max')
    
    # with gridspec:
    import matplotlib.gridspec 
    import matplotlib.pyplot as plt
    fig=plt.figure()

    gs = gridspec.GridSpec(2,2)
    gs.update(right=0.87)
    for ii in range(2):
        for jj in range(2):
            ax=plt.subplot(gs[ii,jj])
            pcm=ax.pcolormesh(rand(10,10))
    colorbarRight(pcm,gs,fig,extend='max')
    '''
    import numpy as np
    import matplotlib.gridspec as gs
    import matplotlib.pyplot as plt
    if type(ax) is gs.GridSpec:
        # gridspecs are different than axes:
        pos = ax.get_grid_positions(fig)
        y0 = pos[0][-1]
        y1 = pos[1][0]
        x1 = pos[3][-1]
    else: 
        if ~(type(ax) is np.ndarray):
            # these are supposedly axes:
            ax=np.array(ax)
        # get max x1, min y0 and max y1
        y1 = 0.
        y0 = 1.
        x1=0.
        for axx in ax.flatten():
            pos=axx.get_position()
            x1=np.max([pos.x1,x1])
            y1=np.max([pos.y1,y1])
            y0=np.min([pos.y0,y0])
    height = y1-y0
    pos2 = [x1 + gap, y0 + (1.-shrink)*height/2.,  width, height*shrink]
    cax=plt.axes(position=pos2)
    fig.colorbar(pcm,cax=cax,**kwargs)

def getCmap2():
    cmap = plt.get_cmap('RdBu_r')
    colors = cmap(np.linspace(0.5, 1, cmap.N // 2))
    # Create a new colormap from those colors
    cmap2 = LinearSegmentedColormap.from_list('Upper Half', colors)
    return cmap2

def plotFluxes(ax,D,dx=50,dy=50):

    H=D['Depth'];x=D['x']/1e3;y=D['y']/1e3
   
    cmap2=getCmap2()
    F=np.abs(D['uPbc']+1j*D['vPbc']); H = D['Depth']; Fu=D['uPbc'];Fv=D['vPbc']
    pcm1=ax.pcolormesh(x,y,F,cmap=cmap2,rasterized=True)
    pcm1.set_clim(0,3.)
    
    xg = np.arange(min(x),max(x),dx);yg = np.arange(min(y),max(y),dy)
    X,Y=np.meshgrid(xg,yg)
    # get Pug and Pvg
    Pug=0.*X;Pvg=0.*Y
    for j in range(np.size(yg)):
        indy = np.where(y>yg[j])[0][0]
        aa=np.interp(xg,x,Fu[indy,:])
        Pug[j,:]=aa
        Pvg[j,:]=np.interp(xg,x,Fv[indy,:])
    ax.quiver(xg,yg,Pug,Pvg,scale=20.,color='0.5',edgecolors='0.5',linewidths=(1,))
    plotdepthcont(ax,x,y,H)
    ax.set_aspect(1.)
    return pcm1
