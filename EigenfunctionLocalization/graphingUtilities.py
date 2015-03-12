import pylab
import matplotlib.pyplot as plt
from fractions import gcd
import numpy as np
import utilities as utes
import itertools as it

if utes.isWin:#Hack for importing mayavi only when on my computer
              #instead of remote linux server
    from mayavi import mlab

import matplotlib.animation as animation

#take a dict of functions and superimpose them
def superImpose(g):
    K = g.keys()
    gnew = np.zeros_like(g[K[0]])
    n = gnew.shape[0]
    for k in K:
        for (i,j) in it.product(xrange(n), repeat = 2):
            if g[k][i,j] != 0:
                gnew[i,j] = g[k][i,j]
    return gnew

#plot a function in 2D
def plotFunction(f, U = None):
    if type(f) == dict:
        f = superImpose(f)
    else:
        f = f.copy()
    M = max(abs(f.flatten())) * (-1)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlim(left = 0, right = 1)
    n = f.shape[0]
    if U == None:
        for i in [0, n-1]:
            for j in xrange(n):
                f[i,j] = M
                f[j,i] = M
    else:
        for i,j in it.product(range(n), repeat = 2):
            if (i, j) not in U.Domain:
                f[i,j] = M
    ax.imshow(f.transpose(),cmap=plt.cm.spectral,origin='lower',extent=(0,1,0,1))
    plt.show()

#plot a function in 3D
def plotFunction3D(f):
    if type(f) == dict:
        f = superImpose(f)
    n = max(f.shape)
    meshSize = 1./n
    x, y = np.mgrid[0:n*meshSize:meshSize, 0:n*meshSize:meshSize]    
    s = mlab.surf(x, y, f, warp_scale = "auto")
    #mlab.axes()
    #cs = contour_surf(x, y, f, contour_z=0)
    return s
    
#for quantumSquare module, plot a potential function in 3D
def plotPotential3D(f):
    M = max(f.flatten())
    n = np.shape(f)[0]
    scale = float(n)/(1.5*M)
    x, y = np.mgrid[0:n:1, 0:n:1]    
    s = mlab.barchart(x, y, f*scale)
    mlab.axes()
    return s
    
def returnPartition(U, maximum = None):
    n = U.n
    h = U.meshSize
    if maximum == None:
        maximum = max(U.keys())
    M = (maximum + 1) * 1.1
    f = np.zeros((n,n)) + M
    for s in U.bDomain:
        f[s] = 0
    for s in U.Domain:
        f[s] = 0
    for key, S in U.Omega.iteritems():
        for s in S:
            f[s] = key + 1
    X = np.array(range(n))
    Y = np.array(range(n))
    X = X*h
    Y = Y*h
    X,Y = np.meshgrid(X,Y)
    f = f.transpose()
    return pylab.imshow(f, cmap = pylab.cm.spectral, origin='lower', extent = (0,1,0,1))
    
#plot a partition; black = region not incorporated
def showPartition(U):
    returnPartition(U)
    pylab.show()
    
#show a partition with mlab;
################################
#WARNING: Not thoroughly tested#
################################
def showPartitionMlab(U):
    n = U.n
    h = U.meshSize
    f = np.zeros((n,n))
    for key, S in U.Omega.iteritems():
        for s in S:
            f[s] = key + 1
    X = np.array(range(n))
    Y = np.array(range(n))
    X = X*h
    Y = Y*h
    X,Y = np.meshgrid(X,Y)
    f = f.transpose()
    mlab.imshow(f, colormap = 'spectral')
    mlab.show()

def showPartOverFxn(U, f):
    if type(f) == dict:
        f = superImpose(f)
    M = max(abs(f.flatten())) * (-1.1)
    fig = plt.figure()
    ax = fig.add_subplot(121)
    ax.set_xlim(left = 0, right = 1)
    n = U.n
    h = U.meshSize
    g = np.zeros((n,n))
    for key, S in U.Omega.iteritems():
        for s in S:
            g[s] = key + 1
    ax.imshow(g.transpose(), cmap=plt.cm.spectral, origin='lower', extent=(0,1,0,1))
    ax = fig.add_subplot(122)
    glines = np.zeros((n,n))
    for key, S in U.bOmega.iteritems():
        for s in S:
            glines[s] = M
    nf = f.shape[0]
    m = gcd(n - 2, nf - 2)
    L = ((n-2)*(nf - 2)) / m
    glines = np.delete(glines, [0, n-1], 0)
    glines = np.delete(glines, [0, n-1], 1)
    glines = utes.expand(glines, (nf-2)/m)
    flines = np.delete(f, [0, nf - 1], 0)
    flines = np.delete(flines, [0, nf-1], 1)
    flines = utes.expand(flines, (n - 2)/m)
    for I in xrange((L + 2)**2):
        i, j = I/(L+2), I % (L+2)
        if glines[i,j] != 0:
            flines[i,j] = M
    #print glines.shape, flines.shape, m, L, n, nf
    ax.imshow(flines.transpose(),cmap=plt.cm.spectral,origin='lower',extent=(0,1,0,1))
    plt.show()
     
def showPartOverFxnList(U, flist, xran, yran, scale=False):
    fig = plt.figure()
    I = len(flist) + 1
    ax = fig.add_subplot(xran,yran,1)
    ax.set_xlim(left = 0, right = 1)
    n = U.n
    #h = U.meshsize
    g = np.zeros((n,n))
    for key, S in U.Omega.iteritems():
        for s in S:
            g[s] = key + 1
    ax.imshow(g.transpose(), cmap=plt.cm.spectral, origin='lower', extent=(0,1,0,1))
    for i in xrange(xran*yran-1):      
        f = flist[i]     
        M = max(abs(f.flatten())) * (-1.1)
        #print i + 2
        ax = fig.add_subplot(xran,yran,i+2)
        #print i + 2
        glines = np.zeros((n,n))
        for key, S in U.bOmega.iteritems():
            for s in S:
                glines[s] = M
        nf = f.shape[0]
        m = gcd(n - 2, nf - 2)
        L = ((n-2)*(nf - 2)) / m
        glines = np.delete(glines, [0, n-1], 0)
        glines = np.delete(glines, [0, n-1], 1)
        glines = utes.expand(glines, (nf-2)/m)
        flines = np.delete(f, [0, nf - 1], 0)
        flines = np.delete(flines, [0, nf-1], 1)
        flines = utes.expand(flines, (n - 2)/m)
        for I in xrange((L + 2)**2):
            i, j = I/(L+2), I % (L+2)
            if glines[i,j] != 0:
                flines[i,j] = M
        #print glines.shape, flines.shape, m, L, n, nf
        #kwargs = {'XTick':.5, 'YTick':.5}
        ax.imshow(flines.transpose(),cmap=plt.cm.spectral,origin='lower',extent=(0,1,0,1))
    plt.show()
    
def makeMovie(Wlist, name, fps = 15):
    """This is system dependent; must have ffmpeg installed."""
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps = fps, metadata = dict(artist='Me'), bitrate = 1800)
    fig1 = plt.figure()
    #plt.xlim(0, 1)
    #plt.ylim(0, 1)
    #plt.xlabel('x')
    #plt.title('evolution')
    maximum = -1
    for W in Wlist:
        maximum = max(maximum, W.maxKey())
    frames = [(returnPartition(w, maximum = maximum),) for w in Wlist]
    im_ani = animation.ArtistAnimation(fig1, frames, interval=50, repeat_delay=3000,
        blit=True)
    im_ani.save(name+'.mp4', writer=writer)