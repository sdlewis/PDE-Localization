import math
import numpy as np
from twoDPartition import partition
from scipy.sparse.linalg import spsolve, eigsh, minres, cg
from scipy.linalg import eigh as sceigh
from scipy.linalg import solve
import utilities as utes
import graphingUtilities as gUtes
import setup
import time
from mpl_toolkits.mplot3d import Axes3D


#For testing purposes
#import time

#Global parameters set here

#creating a matrix to represent an operator
opToMatrix = setup.operatorToMatrix 

#Initial temperature used in simulated annealing
Tinit = setup.Tinit

#The linear and square term on the volume functional
Clinear = setup.Clinear
Csquare = setup.Csquare

#Sparse Solvers
tol = pow(10, -7)
sparseSolve = spsolve

#The volume cost function
def G(x): return Clinear * x + Csquare * pow(x,2)

#The energy due to volume
def volumeEnergy(U, index = None):
    if index == None:
        return volumeEnergy(U, U.keys())
    elif type(index) == list or type(index) == set:
        return {i: volumeEnergy(U, i) for i in index}
    elif type(index) == int:
        return G(U.volume(index))
    else:
        raise ValueError, "volumeEnergy index takes None, list/set, or int"

#the total energy
def energy(U, f, index = None):
    if index == None:
        return energy(U, f, U.keys())
    elif type(index) == list or type(index) == set:
        return {i: energy(U, f, i) for i in index}
    elif type(index) == int:
        if type(f) == dict:
            return energy(U, f[index], index)
        elif type(f) == np.ndarray:
            return setup.functionEnergy(U, f, index) + volumeEnergy(U, index)
        else:
            raise ValueError, "f must be np.array or dict of np.arrays in energy"
    else:
        raise ValueError, "volumeEnergy index takes None, list/set, \'domain\',or int"     
        
        
###############################################
#                                             #
# Integrate code for the landscape solving    #
#                                             #
###############################################

#Solve for the landscape function on the component U.Omega[i]

def landscape(U, index = None):
    if index == None:
        return landscape(U, U.keys())
    elif type(index) == list or type(index) == set:
        return {i:landscape(U, i) for i in index}
    elif type(index) == int or index == 'domain':
        if index in U.keys() or index == 'domain':
            M = U.localOperator(index)
            k = M.shape[0]
            q = np.zeros(k) + U.meshSize**2
            if k == 1:
                u = solve(M.todense(), q)
            else:
                #u, flag = sparseSolve(M, q, tol=tol)
                u = sparseSolve(M, q)
            u = utes.unindex(U, u, index)
            #if flag == 1:
            #    print "Error in landscape linear algebra"
            return u
        else:
            return np.zeros((U.n, U.n))
    else:
        raise ValueError, "index in landscape must be None, list, set or int"

#Solve for the landscape function on the component U.Omega[i]
#archived as of 3/19/14
def landscapeOLD(U, index = None):
    if index == None:
        return landscapeOLD(U, U.keys())
    elif type(index) == list or type(index) == set:
        return {i:landscapeOLD(U, i) for i in index}
    elif type(index) == int or index == 'domain':
        if index in U.keys() or index == 'domain':
            M, indexing = opToMatrix(U, index)
            q = np.zeros(len(indexing)) + U.meshSize**2
            if len(indexing) == 1:
                u = solve(M.todense(), q)
            else:
                u = spsolve(M.tocsr(), q)
            u = utes.unindexer(U, u, indexing)
            return u
        else:
            return np.zeros((U.n, U.n))
    else:
        raise ValueError, "index in landscapeSolver must be None, list, set or int"


#Find the landscape function over the domain        
def globalLandscape(U):
    W = U.copy()
    W.Omega = {0: W.Domain}
    W.updateBoundary()
    return landscape(W, 0)
    
def globalLandscape2(U):
    return landscape(U, 'domain')
    
#Returns true iff deltaJ is negative (i.e., only when the functional improves)
def lessthan(deltaJ, total_time, iterations):
    return deltaJ < 0

#simulated annealing decision method
def simulatedAnnealing(deltaJ, time, iterations, Tinit = Tinit):
    T = (Tinit *(iterations - time))/iterations
    return P(deltaJ, T) > np.random.rand()
    
#function for simulatedAnnealing method
#Standard function, though others exist
def P(deltah, T):
    if deltah < 0:
        return 1
    else:
        return math.exp((-deltah)/T)

#perform a search for a minimum of the functional
#fps = frames per second; ipf = iterations per frame
def runSearch(U, indices = None, decisionMethod = lessthan, 
                iterations = setup.iterations, makeMovie = False, ipf = 20,
                initMerge = 150, mergePeriod = 200, printUpdate = 50):
                
    tlast = time.time()
    t0 = tlast
                
    UCurr = U.copy()
    UBest = UCurr.copy()
    
    if makeMovie:
        UMovieList = list()
                
    if indices == None:
        indices = U.keys()
    elif type(indices) != set:
        raise ValueError, "indices must be a set in runSearch"
        
    uCurr = landscape(UCurr, indices)
    JCurr = energy(UCurr, uCurr)
    JBest = JCurr.copy()
    totalTime = 0
    
    tlast = time.time()-tlast
    print "setup time = %s secs" % tlast
    
    while totalTime < iterations:
        for index in indices:
            if index in UCurr.keys():
                UNew = UCurr.copy()
          
                #print 'UNew', UNew.keys()
                uNew = uCurr.copy()
                updates = UNew.evolve(index, returnUpdates = True)
                #print 'UNew 2', UNew.keys()
                legalUpdates = updates.intersection(UNew.keys())
                uUpdates = landscape(UNew, legalUpdates)

                (uNew.pop(i) for i in uNew if i not in UNew)
                uNew.update(uUpdates)      
                
                JNew = energy(UNew, uNew)
                deltaJ = sum(JNew.get(i,0) - JCurr.get(i,0) for i in JCurr.iterkeys())
                
                if decisionMethod(deltaJ, totalTime, iterations):
                    #print 'decisionMade', sum(JCurr.values()), sum(JNew.values())
                    #print JCurr, JNew
                    uCurr.update(uUpdates)
                    UCurr = UNew
                    JCurr = JNew
                    if sum(JCurr.values()) < sum(JBest.values()):
                        JBest = JCurr.copy()
                        UBest = UCurr.copy()
                        uBest = uCurr.copy()
                        
        totalTime += 1
        if makeMovie and totalTime%ipf == (1%ipf):
            UMovieList.append(UCurr)
            
        if totalTime % printUpdate == 1:
            tlast = int(time.time()-t0)
            mins = tlast/60
            secs = tlast %60
            data = (sum(JCurr.values()), sum(JBest.values()), totalTime, 
                    100.*totalTime/iterations, mins, secs)
            print 'curr=%s, best=%s,\ntotalTime=%s, %s/100, %s mins, %s secs'%data
        
        #Merging routine    
        #TO DO: Test 
        if ((totalTime % mergePeriod == (mergePeriod/2 + initMerge)%mergePeriod
            and totalTime > initMerge) or (totalTime == initMerge)):   
            isImproving = True
            print 'merging...'
            while isImproving:
                deltas = dict()
                relDeltas = dict()
                JMerge = dict()
                uMerge = dict()
                edges = UCurr.adjGraph()
                #print 'edges are ', edges.edges
                if len(edges) == 0:
                    break                
                for e in edges:
                    if e not in deltas:
                        UMerge = UCurr.copy()
                        UMerge.merge(e[0], e[1])
                        uMerge[e] = landscape(UMerge, index = e[0])
                        JMerge[e] = energy(UMerge, uMerge[e], e[0])
                        deltas[e] = JMerge[e] - JCurr[e[0]] - JCurr[e[1]]
                        relDeltas[e] = deltas[e] /abs(JCurr[e[0]] + JCurr[e[1]])
                        #print 'loop', JMerge, deltas
                minKey = min(deltas, key = lambda x: deltas[x])
                relativeMinKey = min(relDeltas, key = lambda x: deltas[x])
                #print minKey, deltas[minKey]
                mergeKey = relativeMinKey
                if deltas[minKey] < 0:
                    UCurr.merge(mergeKey[0], mergeKey[1])
                    JCurr.pop(mergeKey[1])
                    JCurr[mergeKey[0]] = JMerge[mergeKey]
                    uCurr[mergeKey[0]] = uMerge[mergeKey]
                    K = deltas.keys()
                    #print 'keys are ', K
                    #print 'merging ', minKey
                    for e in K:
                        #print 'e is ', e
                        if mergeKey[0] in e or mergeKey[1] in e:
                            deltas.pop(e)
                            #print 'popping ', e
                else:
                    isImproving = False
                                
    
    returnValues = {'best': UBest, 'curr':UCurr}
    
    if makeMovie:
        UMovieList.append(UBest)
        returnValues['movie'] = UMovieList
    
    return returnValues
    

    
if __name__ == '__main__':
    n = setup.n
    U = partition(n, 1./n)
    U.seed(10)
    
    
    
    #S = list()
    #S.append(utes.graphball((80,80), 40, U.Domain))
    #S.append(utes.graphball((40,20), 9, U.Domain))
    #S.append(utes.graphball((20,40), 9, U.Domain))
    #S.append(utes.graphball((40,40), 9, U.Domain))
    #for s in S:
    #    U.add(s)
    g = landscape(U)
    UCurr = U.copy()
    uCurr = landscape(UCurr)
    JCurr = energy(UCurr, uCurr)
    print 'JCurr', JCurr

    