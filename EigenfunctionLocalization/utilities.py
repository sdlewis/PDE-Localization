import numpy as np
import pickle
import math
import os
import inspect

if os.sys.platform[:3].lower() == 'win':
    isWin = True
    slash = '\\'
else:
    isWin = False
    slash = '/'

cwd = inspect.getfile(inspect.currentframe())[:-12]
dataFolder = cwd + '..' + slash + '..' + slash + 'data' + slash
os.chdir(dataFolder)

#Saving and Loading with Pickle
def picklesave(U, name):
    f = open(dataFolder + name, 'wb')
    pickle.dump(U, f)
    f.close()
    
def pickleload(name):
    f = open(dataFolder + name, 'rb')
    U = pickle.load(f)
    f.close()
    return U

#return the neighbors of a point or set of points along axes
def neighbors(points):
    """return the manhattan neighbors of a point stored as a tuple
    or set/list of points stored as tuples
    """
    deltas = ((1,0), (-1,0), (0,1), (0,-1))
    if type(points) == tuple:
        x1, x2 = points
        return set([(x1+d1, x2+d2) for (d1, d2) in deltas])
    elif type(points) == list or type(points) == set:
        S = set()
        for d1, d2 in deltas:
            for x1, x2 in  points:
                S.add((x1+d1, x2+d2))
        return S
    else:
        raise ValueError, "tuple(point), list, or set only in neighbors"
    
#Provide linear indexing of interior points
def indexer(U, index):
    if index == 'domain':
        m = 0
        indexing = dict()
        for (i,j) in U.Domain:
            indexing[(i,j)] = m
            m += 1
        return indexing      
    elif index not in U.keys():
        print U.keys()
        print U.Omega
        raise ValueError, 'index' + str(index)+  ' not in U.keys()'
        
    m = 0
    indexing = dict()
    for (i,j) in U.Omega[index]:
        indexing[(i,j)] = m
        m += 1
    return indexing
    
#unindex a function to be defined on the nxn square
def unindex(U, u, index):
    n = U.n
    u_ = np.zeros((n,n))
    for (a,v) in zip(U[index], u):
        u_[a] = v
    return u_
    
#unindex a function to be defined on the nxn square
#Deprecated 3/19/14
def unindexer(U, u, indexing):
    n = U.n
    u_ = np.zeros((n,n))
    for key, value in indexing.iteritems():
        u_[key] = u[value]
    return u_
    
#Fast Matrix Slicing
def cooSlice(M, indices):
    M2 = M.tocsr()
    M2 = M2[indices]
    #M2 = M2.tocsc() -- Optional:
    #Doesn't seem to affect performance
    M2 = M2[:,indices]
    return M2
    
def csrSlice(M,indices):
    return M[indices][:,indices]
    
def csrSlice2(M, indices):
    M2 = M[indices]
    M2 = M2.tocsc()
    return M2[:, indices]
    
#return the neighbors of a point or set of points along axes and diagonals
def extendedNeighbors(points):
    deltas = ((1,0), (-1,0), (0,1), (0,-1), (1,1), (-1,-1), (1,-1), (-1,1))
    if type(points) == tuple:
        x1, x2 = points
        return set([(x1+d1, x2+d2) for (d1, d2) in deltas])
    elif type(points) == list or type(points) == set:
        S = set()
        for d1, d2 in deltas:
            for x1, x2 in  points:
                S.add((x1+d1, x2+d2))
        return S
    else:
        raise ValueError, "tuple(point), list, or set only in extendNeighbors"
   
#Find connected component in extended neighbor graph
def extendedNeighborComponent(v, V):    
    component = set([v])
    N = set([v])
    while len(N) > 0:
        N = extendedNeighbors(N).intersection(V)
        component.update(N)
        V.difference_update(N)
    return component 
   
#return the points r close to N in S by utes.extendedNeighbors' graph metric
def graphball(N, r, S):
    if type(N) == tuple:
        N = set([N])
    elif type(N) == set:
        N = N.copy()
    else:
        raise ValueError, "graphball takes N as set, or tuple (point)"
    V = S.copy()
    B = N.copy()
    while len(N) > 0 and r > 0:
        N = extendedNeighbors(N).intersection(V)
        B.update(N)
        V.difference_update(N)
        r -= 1
    return B
    

def euclideanBall(N, r, S):
    if type(N) == tuple:
        N = set([N])
    elif type(N) == set:
        N = N.copy()
    else:
        raise ValueError, "graphball takes N as set, or tuple (point)"
    B = set()
    for s in S:
        for n in N:
            if sum((x-y)**2 for x,y in zip(s, n))**.5 < r:
                B.add(s)
                break
    return B
    
#Takes function V as a kxk input
#expands each cell by l and adds boundary values
#end size is nxn for n = l*k + 2
def expand(f, l):
    k = np.shape(f)[0]
    n = l*k + 2
    V = np.zeros((n,n))
    for I in xrange(k**2):
        i, j = I/k, I%k
        for L in xrange(l**2):
            di, dj = L/l + 1, L%l +1
            V[i*l + di, j*l + dj] = f[i,j]
    return V
    
#energy methods
#def gradSquare1(u, U, index):
#    """Return the integral of the gradient squared over U[i]"""
#    total = 0
#    for (i,j) in U.Omega[index]:
#        total += pow(u[i,j] - u[i+1, j],2) + pow(u[i,j] - u[i, j+1],2)
#    return total

def gradSquare(u,U=None,index=None):
    ucent = u[:-1, :-1]
    uleft = u[1:, :-1]
    udown = u[:-1, 1:]
    dx = (ucent - uleft)**2
    dy = (ucent - udown)**2
    return dx.sum() + dy.sum()

def integrate1(u, U, index):
    return U.meshSize**2 * sum(u[p] for p in U.Omega[index])

def integrate(u, U, index=None):
    return (u.sum())*(U.meshSize**2)
    
#return a uniform random point in a set
def randomPoint(S):
    """Return a uniform random point in a set S"""
    Slist = list(S)
    l = len(Slist)
    i = int(np.random.rand()*l)
    return Slist[i] 
