import numpy as np
from scipy.sparse import coo_matrix
import utilities as utes
    
#Declared Variables

n = 75
  
iterations = 201
Tinit = pow(10, -11)
Clinear = 0.000001
Csquare = .0035

def laplaceToMatrix(U, index):

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
    data = np.array([4 for p in indexing.iterkeys()])
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
    
def laplaceToMatrixOLD(U, index):
    """Laplacian matrix constructor"""
    
    indexing = utes.indexer(U, index)
    l = len(U.Omega[index])
    row = np.array([indexing[p] for p in indexing.iterkeys()])
    col = np.array([indexing[p] for p in indexing.iterkeys()])
    data = np.array([4 for p in indexing.iterkeys()])
    row2 = np.zeros(4*l)
    col2 = np.zeros(4*l)
    data2 = np.zeros(4*l)
    i = 0
    for p in indexing.iterkeys():
        if p in U.abOmega[index]:
            for q in utes.neighbors(p):
                if q not in U.bOmega[index]:
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
    
operatorToMatrix = laplaceToMatrix

def globalOperator(U):
    return operatorToMatrix(U,'domain')
    
#Energy Methods
def __functionEnergy(U, f, index):
    energy = .5 * utes.gradSquare(f, U, index)
    energy -= utes.integrate(f, U, index)
    return energy
    
def functionEnergy(U, f, index):
    return __functionEnergy(U, f, index)