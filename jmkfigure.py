#!/usr/local/bin/python
# Filename: jmkfigure.py

# from matplotlib import rc
from pylab import *

def djmkfigure(width,vext):
    """
    djmkfigure(width,vext):
    width is column widths, and vext is fractional 10 page height.  
    """
    wid = 3*width+3./8.;
    height = 10*vext;
    rc('figure',figsize=(wid,height),dpi=96)
    rc('font',size=9)
    rc('font',family='sans-serif');
    # rcParams['font.sans-serif'] = ['Verdana']
    rc('axes',labelsize='large') 
    leftin = 0.75
    rightin = 0.25
    botin = 0.4
    rc('figure.subplot',left=leftin/wid) 
    rc('figure.subplot',right=(1-rightin/wid)) 
    rc('figure.subplot',bottom=botin/height) 

def jmkprint(fname,pyname,dirname='doc'):
    """
    def jmkprint(fname,pyname)
    def jmkprint(fname,pyname,dirname='doc')
    """
    import os
    
    try:
        os.mkdir(dirname)
    except:
        pass

    if dirname=='doc':
        pwd=os.getcwd()+'/doc/'
    else:
        pwd=dirname+'/'
    savefig(dirname+'/'+fname+'.pdf',dpi=400)
    savefig(dirname+'/'+fname+'.png',dpi=400)
    
    fout = open(dirname+'/'+fname+'.tex','w')
    str="""\\begin{{figure*}}[htbp]
  \\begin{{center}}
    \\includegraphics[width=\\twowidth]{{{fname}}}
    \\caption{{
      \\tempS{{\\footnotesize {pwd}/{pyname} ;     
        {pwd}{fname}.pdf}}
      \\label{{fig:{fname}}} }}
  \\end{{center}}
\\end{{figure*}}""".format(pwd=pwd,pyname=pyname,fname=fname)
    fout.write(str)
    fout.close()
    
    cmd = 'less '+dirname+'/%s.tex | pbcopy' % fname
    os.system(cmd) 


def tsdiagramjmk(salt,temp,cls=[]):
    import numpy as np
    import seawater
    import matplotlib.pyplot as plt
     
     
    # Figure out boudaries (mins and maxs)
    smin = salt.min() - (0.01 * salt.min())
    smax = salt.max() + (0.01 * salt.max())
    tmin = temp.min() - (0.1 * temp.max())
    tmax = temp.max() + (0.1 * temp.max())
 
    # Calculate how many gridcells we need in the x and y dimensions
    xdim = round((smax-smin)/0.1+1,0)
    ydim = round((tmax-tmin)+1,0)

     
    # Create empty grid of zeros
    dens = np.zeros((ydim,xdim))
     
    # Create temp and salt vectors of appropiate dimensions
    ti = np.linspace(1,ydim-1,ydim)+tmin
    si = np.linspace(1,xdim-1,xdim)*0.1+smin
     
    # Loop to fill in grid with densities
    for j in range(0,int(ydim)):
        for i in range(0, int(xdim)):
            dens[j,i]=seawater.dens(si[i],ti[j],0)
     
    # Substract 1000 to convert to sigma-t
    dens = dens - 1000
 
    # Plot data ***********************************************
    if not(cls==[]):
        CS = plt.contour(si,ti,dens, cls,linestyles='dashed', colors='k')
    else:
        CS = plt.contour(si,ti,dens,linestyles='dashed', colors='k')
        
    plt.clabel(CS, fontsize=9, inline=1, fmt='%1.2f') # Label every second level
    ax1=gca()
    #    ax1.plot(salt,temp,'or',markersize=4)
     
    ax1.set_xlabel('S [psu]')
    ax1.set_ylabel('T [C]')

#########################################################################

def facetpcolor(x,y,z,**kwargs):
    # keep y the same.  Expland x
    x0=1.*x
    [M,N]=shape(z)
    dx = diff(x)
    x = np.tile(x,(2,1))
    x[0,1:]=x[0,1:]-dx/2.
    x[0,0]=x[0,0]-dx[0]/2.
    x[1,:-1]=x[1,:-1]+dx/2
    x[1,-1]=x[1,-1]+dx[-1]/2

    z = np.tile(z,(1,1,2))
    [m,n]=shape(x)
    x=np.reshape(x.T,[2*n],order='C')
    z=np.reshape(z,[M,-1],order='C')
    zz = 0.*z
    zz[:,0:-1:2]=z[:,0:N]
    zz[:,1:-1:2]=z[:,(N):(2*N-1)]
    
    return pcolormesh(x,y,zz,**kwargs) 

#########################################################################

def pcolormeshRdBu(x,y,z,**kwargs):
    # def pcolormeshRdBu(x,y,z,**kwargs):
    return pcolormesh(x,y,z,rasterized=True,cmap=cm.RdBu_r,**kwargs)    

####################
def gmtColormap(fileName,GMTPath = '~/python/cmaps/'):
      """
      gmtColormap(fileName,GMTPath='~/python/cmaps/')
      
      Returns a dict for use w/ LinearSegmentedColormap

      cdict = gmtcolormapPylab.gmtcolormapPylab('spectrum-light')
      colormap = pylab.cm.colors.LinearSegmentedColormap('spectrum-light',cdict)
      """
      import colorsys
      import numpy as N
      filePath = GMTPath+"/"+ fileName +".cpt"
      try:
          f = open(filePath)
      except:
          print "file ",filePath, "not found"
          return None

      lines = f.readlines()
      f.close()

      x = N.array([])
      r = N.array([])
      g = N.array([])
      b = N.array([])
      colorModel = "RGB"
      for l in lines:
          ls = l.split()
          if l[0] == "#":
             if ls[-1] == "HSV":
                 colorModel = "HSV"
                 continue
             else:
                 continue
          if ls[0] == "B" or ls[0] == "F" or ls[0] == "N":
             pass
          else:
              x=N.append(x,float(ls[0]))
              r=N.append(r,float(ls[1]))
              g=N.append(g,float(ls[2]))
              b=N.append(b,float(ls[3]))
              xtemp = float(ls[4])
              rtemp = float(ls[5])
              gtemp = float(ls[6])
              btemp = float(ls[7])
              

      x=N.append(x,xtemp)
      r=N.append(r,rtemp)
      g=N.append(g,gtemp)
      b=N.append(b,btemp)

      nTable = len(r)
      if colorModel == "HSV":
         for i in range(r.shape[0]):
             rr,gg,bb = colorsys.hsv_to_rgb(r[i]/360.,g[i],b[i])
             r[i] = rr ; g[i] = gg ; b[i] = bb
      if colorModel == "HSV":
         for i in range(r.shape[0]):
             rr,gg,bb = colorsys.hsv_to_rgb(r[i]/360.,g[i],b[i])
             r[i] = rr ; g[i] = gg ; b[i] = bb
      if colorModel == "RGB":
          r = r/255.
          g = g/255.
          b = b/255.
      print shape(x)
      xNorm = (x - x[0])/(x[-1] - x[0])

      red = []
      blue = []
      green = []
      for i in range(len(x)):
          red.append([xNorm[i],r[i],r[i]])
          green.append([xNorm[i],g[i],g[i]])
          blue.append([xNorm[i],b[i],b[i]])
      colorDict = {"red":red, "green":green, "blue":blue}
      return (colorDict)    
    
####################
def gmtcmap(fileName,GMTPath  = '/Users/jklymak/python/cmaps/'):
    import matplotlib.pylab as pylab
    
    cdict=gmtColormap(fileName,GMTPath)
    colormap = pylab.cm.colors.LinearSegmentedColormap(fileName,cdict)
    return colormap
          
#####################
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

    if type(ax) is matplotlib.gridspec.GridSpec:
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
    cax=axes(position=pos2)
    fig.colorbar(pcm,cax=cax,**kwargs)
