import numpy as np
from scipy.sparse import coo_matrix
import utilities as utes

#Takes potential function V as a kxk input
#expands each cell by k and adds boundary values
#end size is nxn for n = l*k + 2
def expandV(V_, l):
    return utes.expand(V_, l)

#Assumes k|(n-2)
#Returns a kxk random grid, each square uniformly
#selected from 0 to 8000. The grid is then fitted
#to an nxn square with 0s on the boundary. 
#returns both the fitted matrix V and the unfitted V_
def randomV(n,k):
    m = (n-2)/k
    V_ = np.random.rand(k,k)*8000
    V = expandV(V_, m)
    return V, V_

#Declared Variables

#k = 3 #number of ticks for the potential functions
#m = 4 #number of meshpoint ticks ****PER**** potential function tick
#n = k*m+2 #Will be U.n in
#V, V_ = randomV(n,k)

V_ = utes.pickleload('V.dat')
k = V_.shape[0]
m = 8
n = k*m + 2
V = expandV(V_, m)

seeds = 400
iterations = 2501
Tinit = pow(10, -12)
Clinear = .000002
Csquare = .00004

#Returns M, matrix of discretized operator and indexing from indexer
def operatorToMatrix(U, index):
    """QuantumSquare matrix constructor"""
    
    indexing = utes.indexer(U, index)
    h = U.meshSize
    hsq = pow(h,2)
    
    if index == 'domain':
        dom = U.Domain
        bdom = U.bDomain
        abdom = U.abDomain
    else:
        dom = U[index]
        bdom = U.bOmega[index]
        abdom = U.abOmega[index]
    
    l = len(dom)
    row = np.array([indexing[p] for p in indexing.iterkeys()])
    col = np.array([indexing[p] for p in indexing.iterkeys()])
    data = np.array([4 + hsq * V[p] for p in indexing.iterkeys()])
    row2 = np.zeros(4*l)
    col2 = np.zeros(4*l)
    data2 = np.zeros(4*l)
    i = 0
    for p in indexing.iterkeys():
        if p in abdom:
            for q in utes.neighbors(p):
                if q not in bdom:
                    row2[i] = indexing[p]
                    col2[i] = indexing[q]
                    data2[i] = - 1.
                    i += 1
        else:
            for q in utes.neighbors(p):
                row2[i] = indexing[p]
                col2[i] = indexing[q]
                data2[i] = - 1.
                i += 1
    data = np.append(data, data2[:i])
    row = np.append(row, row2[:i])
    col = np.append(col, col2[:i])
    M = coo_matrix((data, (row, col)), shape = (l, l))
    return M, indexing
    
def globalOperator(U):
    return operatorToMatrix(U,'domain')

#Energy Methods
def __functionEnergy(U, f, V, index):
    energy = .5 * utes.gradSquare(f, U, index)
    energy += utes.integrate(V*(f**2), U, index)
    energy -= utes.integrate(f, U, index)
    return energy
    
def functionEnergy(U, f, index):
    return __functionEnergy(U, f, V, index)
    

