import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import time
import numpy as np
import cPickle as pickle

img=mpimg.imread('../ttide15/DS96F4.png')

fig,ax=plt.subplots()
ax.imshow(img)
fig.show()
x0v=np.float(raw_input("x lower left"))
x1v=np.float(raw_input("x lower right"))
y0v=np.float(raw_input("y lower left"))
y1v=np.float(raw_input("y upper left"))

print("lower left")


a=plt.ginput(1)
ax.plot(a[0][0],a[0][1],'x')
fig.canvas.draw()
print("lower right")
b=plt.ginput(1)
ax.plot(b[0][0],b[0][1],'x')
fig.canvas.draw()
xx=np.array([a[0][0],b[0][0]])
xy = np.array([a[0][1],b[0][1]])
print xx,xy
ax.plot(xx,xy,'-')
fig.canvas.draw()
print("click upper left")
c=plt.ginput(1)
yx=np.array([a[0][0],c[0][0]])
yy = np.array([a[0][1],c[0][1]])
print yx,yy
ax.plot(yx,yy,'-')


fig.canvas.draw()

allDone=False
while not(allDone):
    x0 = np.array([])
    y0 = np.array([])
    done = False
    while not(done):
        print('click')
        xxx = plt.ginput(1)
        print(xxx)
        print len(xxx)
        if not(len(xxx)==0):
            print(not(len(xxx)==0))
            x0=np.hstack((x0,xxx[0][0]))        
            y0=np.hstack((y0,xxx[0][1])  )      
            ax.plot(x0,y0,'r')
            fig.canvas.draw()
        else:
            done=True
            name = raw_input('Filename: ')
            lin=dict()
            lin['x0']=x0
            lin['y0']=y0
            ### THIS ASSUMES IMAGE IS SQUARE...
            lin['x'] = (x0-xx[0])/(xx[1]-xx[0])*(x1v-x0v)+x0v
            lin['y'] = (y0-yy[0])/(yy[1]-yy[0])*(y1v-y0v)+y0v
            print lin['y'],lin['x']
            pickle.dump(lin,open(name,'wb'))
            more = raw_input('More lines? [y/n] ')
            if more=='n':
                allDone=True
            else: 
                allDone=False
            
print x0,y0
