import utilities as utes
import graphingUtilities as gUtes
import numpy as np
import itertools as it
from graph import graph, weightedGraph
import setup

#Central object of consideration: partitions

class partition(object):
    def __init__(self, n, meshSize = None, Domain = None, bDomain = None,
                  abDomain = None,
                  Omega = None, bOmega = None, abOmega = None, 
                  runUpdate = False, setOperator = True):
        """n = number of ticks in each dimension, Domain = underlying domain,
        bDomain = boundary of domain, Omega = Partition as a dict,
        bOmega = partition boundaries as a dict, abOmega = adjacent boundary
        runUpdate; if True, bOmega and abOmega will update"""
        
        self.n = n
        if meshSize == None:
            self.meshSize = 1./n
        else:
            self.meshSize = meshSize
        
        #set Domain
        if Domain == None:
            self.Domain = set([(I/(n-2) + 1, I%(n-2) + 1) for I in xrange((n-2)**2)])
        else:
            self.Domain = set(Domain)
            
        #set bDomain
        if bDomain == None:
            self.bDomain = set([(k, 0) for k in xrange(1,n-1)])
            self.bDomain.update([(0,k) for k in xrange(1,n-1)])
            self.bDomain.update([(k, n-1) for k in xrange(1,n-1)])
            self.bDomain.update([(n-1, k) for k in xrange(1,n-1)])
        else:
            self.bDomain = set(bDomain)
            
        #set abDomain
        if abDomain == None:
            self.abDomain = set([(k, 1) for k in xrange(1,n-1)])
            self.abDomain.update([(1,k) for k in xrange(1, n-1)])
            self.abDomain.update([(k, n-2) for k in xrange(1, n-1)])
            self.abDomain.update([(n-2, k) for k in xrange(1, n-1)])
        else:
            self.abDomain = set(abDomain)
            
        #set Omega
        if Omega == None:
            self.Omega = dict()
        else:
            self.Omega = {key: value.copy() for key, value in Omega.iteritems()}
            
        #set bOmega
        if bOmega == None:
            self.bOmega = dict()
        else:
            self.bOmega =  {key: value.copy() for key, value in bOmega.iteritems()}
            
        #set abOmega
        if abOmega == None:
            self.abOmega = dict()
        else:
            self.abOmega =  {key: value.copy() for key, value in abOmega.iteritems()}
            
        #Update the boundary from the sets given
        if runUpdate == True:
            self.updateBoundary()
            
        #set the archive of keys to include the current keys
        self.keyArchive = self.keys().copy()
        
        if setOperator:
            #set the operator and indexing of the domain
            partition.cooOp, partition.indexing = self.globalOperator()
            partition.revIndexing = {v:k for (k,v) in partition.indexing.iteritems()}
            partition.csrOp = partition.cooOp.tocsr()
        
            
    def __iter__(self):
        for i in self.Omega.iterkeys():
            yield i
            
    def __len__(self):
        return len(self.Omega)
        
    def len(self):
        '''Returns the number of keys in self.Omega'''
        return len(self.Omega)
       
    def length(self):
        return self.len()
        
    def __getitem__(self, i):
        if type(i) == tuple:
            key = i[1]
            i = i[0]
        elif i == 'domain':
            return self.Domain
        else:
            key = 'omega'
       # print 'i=%s, key = %s' %(i, key)
        if key[0] == 'o':#omega
            return self.Omega[i]
        if key[0] == 'b':#boundary
            return self.bOmega[i]
        if key[0] == 'c':#omega and its boundary
            return self.completeOmega(i)
            
        print 'nothing to return in partition.__getitem__'\
        
    #########################
    #                       #
    # Simple utilities      #
    #                       #
    #########################
                  
    
    def keys(self):
        return set(self.Omega.keys())
            
    def copy(self):
        '''Create a deep copy of the partition'''
        W = partition(self.n, self.meshSize, self.Domain, self.bDomain, 
            self.abDomain, self.Omega, self.bOmega, self.abOmega,
            setOperator=False)
            
        return W
        
    def maxKey(self):
        '''Return the highest key in the dictionary Omega'''
        if len(self) > 0:
            return max(self.keys())
        else:
            return -1
            
    def newIndex(self):
        '''Generate the minimum new key which has not
        yet been used. Returns 0 if self.keyArchive is 
        empty.'''
        if len(self.keyArchive) == 0:
            return 0
        return max(self.keyArchive) + 1
        
    #########################
    #                       #
    # Geometric data        #
    # retrieval methods     #
    #                       #
    #########################
        
    def completeOmega(self, index):
        '''Return the interior of the ith component of Omega
        together with its boundary mesh points'''
        return self.Omega[index].union(self.bOmega[index])
            
    def adjGraph(self):
        """adjacency graph"""
        G = graph(V = self.keys())
        K = self.keys().copy()
        for i in self.keys():
            K.remove(i)
            for j in K:
                if self.bOmega[i].intersection(self.bOmega[j]):
                    G.add((i,j))
        return G
        
    def wadjGraph(self):
        """Weighted adjacency graph"""
        G = weightedGraph(graph(V = self.keys()))
        K = self.keys().copy()
        for i in self.keys():
            K.remove(i)
            for j in K:
                S = self.bOmega[i].intersection(self.bOmega[j])
                if S:
                    Sratio = 2. * len(S) / (len(self.bOmega[i])+len(self.bOmega[j]))
                    G.add((i,j), Sratio)
        return G
        
    def volume(self, index = None):
        '''Returns the volume at index. If index is not
        specified, returns a complete list of the the volumes'''
        if index == None:
            return [self.volume(i) for i in self.keys()]
        elif type(index) == int:
            return pow(self.meshSize,2)*(len(self.Omega[index])+len(self.bOmega[index])/2.0)
        elif type(index) == list:
            return [self.volume(i) for i in index]
        else:
            raise ValueError, "partition.volume takes Non, int or list"
    def bOmegaComponents(self, index):
        '''Returns a list of the connected components of the boundary
        bOmega[index]. Index must be an integer'''
        verts = self.bOmega[index].copy()
        P = []
        while len(verts)>0:
            v = verts.pop()
            P.append(utes.extendedNeighborComponent(v, verts))
        return P
        
    def abOmegaComponents(self, index):
        '''Returns a list of the connected components of the adjacenct
        boundary points abOmega[index]. Index must be an integer'''
        verts = self.abOmega[index].copy()
        P = []
        while len(verts)>0:
            v = verts.pop()
            P.append(utes.extendedNeighborComponent(v, verts))
        return P
               
    #########################
    #                       #
    # Geometric data        #
    # manipulation methods  #
    #                       #
    #########################
    
    def updateBoundary(self, index = None, ofDomain = False):
        '''Updates the boundary of a region after the underlying set has been
        updated. If index is an integer i, updates the ith index.
        If index is None, updates all sets in Omega. If
        index is a list or a set, update all components in index.
        If ofDomain = True, updates the boundary of the underlying
        Domain'''
        if ofDomain:
            bD = set()
            abD = set()
            for x in self.Domain:
                for y in utes.neighbors(x):
                    if y not in self.Domain:
                        bD.add(y)
                        abD.add(x)
            self.bDomain = bD
            self.abD = abD
        else:
            if index == None:
                self.updateBoundary(self.keys())
            elif type(index) == list or type(index) == set:
                for i in index:
                    self.updateBoundary(i)
            elif type(index) == int:
                self.Omega[index].difference_update(self.bDomain)
                bW = set()
                abW = set()
                #if len(self.Omega[index]) <=5:
                #    print 'in update', self.Omega[index]
                for x in self.Omega[index]:
                    #print 'next level, ', x
                    for y in utes.neighbors(x):
                        if y not in self.Omega[index]:
                            bW.add(y)
                            abW.add(x)
                self.bOmega[index] = bW
                self.abOmega[index] = abW
            else:
                raise ValueError, "index must be list, set, None or int"
        
    def removeKey(self, index):
        '''Removes the key index from the partition'''
        self.Omega.pop(index, None)
        self.bOmega.pop(index, None)
        self.abOmega.pop(index, None)
        
    def grow(self, index, p, r, delta = 1, returnUpdates = False):
        '''Grows self.Omega[index] by adding in all points which are within
        delta of the ball of center p and radius r interesected with
        the boundary'''
        S = utes.graphball(p, r, self.bOmega[index])
        updated = set()
        for i in xrange(delta - 1):
            S = utes.neighbors(S)
            S.intersection_update(self.Domain)
        self.Omega[index].update(S)
        self.updateBoundary(index)
        updated.add(index)
        for i in self.keys():
            if i != index:
                overlap = self.completeOmega(index).intersection(self.Omega[i])
                if len(overlap) > 0:
                    updated.add(i)
                    self.Omega[i].difference_update(overlap)
                    if len(self.Omega[i]) == 0:
                        self.removeKey(i)
                    else:
                        self.updateBoundary(i)
        if returnUpdates == True:
            return updated
        else:
            return None
            
    def shrink(self, index, p, r, delta = 1, returnUpdates = False):
        '''Shrinks self.Omega[index] by removing all points which are within
        delta of the ball of center p and radius r interesected with
        the adjacency boundary'''
        S = utes.graphball(p, r, self.abOmega[index])
        updated = set()
        for i in xrange(delta - 1):
            S = utes.neighbors(S)
            S.intersection_update(self.Domain)
        self.Omega[index].difference_update(S)
        if len(self.Omega[index]) == 0:
            self.removeKey(index)
        else:
            #update boundary as normal
            self.updateBoundary(index)
        updated.add(index)
        if returnUpdates == True:
            return updated
        else:
            return None
        
    def evolve(self, index, pgrow = .5, delta = 1, returnUpdates = False):
        '''Evolve a partition by growing or shring self.Omega[index]
        according to the grow and shrink methods. This occurs with 
        probability of growing == pgrow and a thickness parameter of delta'''
        if np.random.rand() <= pgrow:
            comp = self.bOmegaComponents(index)
            c = len(comp)
            i = int(np.random.rand() * c)
            l = len(comp[i])/2
            r = int(np.random.rand() * l)
            p = utes.randomPoint(comp[i])
            return self.grow(index, p, r, delta = delta, returnUpdates = returnUpdates)
        else:
            comp = self.abOmegaComponents(index)
            c = len(comp)
            i = int(np.random.rand() * c)
            l = len(comp[i])/2
            r = int(np.random.rand() * l)
            p = utes.randomPoint(comp[i])
            return self.shrink(index, p, r, delta = delta, returnUpdates = returnUpdates)
    

    def update(self, W):
        """Update a partition with the data of another partition.
        All keys that are in W but not U will be added. All shared
        keys of U and W will be overwritten.
        
        WARNING: Does NOT check for overlapping sets or incompatabilities.
        Must use updateBoundary separately."""
        if type(W) == partition:
            W = W.copy()
            self.Omega.update(W.Omega)
            self.bOmega.update(W.bOmega)
            self.abOmega.update(W.abOmega)
        else:
            raise ValueError, "partition.update takes W as parition"
        
    def merge(self, target, mergeList):
        """Merge two parts of self together.
        Manages boundaries efficiently. target must be an int
        in self.keys(). mergeList may be an index or a list of indices"""
        if type(mergeList) not in [set, list, tuple]:
            mergeList = set([mergeList])
        if type(mergeList) != set:
            mergeList = set(mergeList)
        newSet = set()
        mergeList.add(target)
        mergeListCopy = mergeList.copy()
        for i in mergeList:
            newSet.update(self.Omega[i])
        for i in mergeList:
            mergeListCopy.remove(i)
            for j in mergeListCopy:
                newSet.update(self.bOmega[i].intersection(self.bOmega[j]))
            self.removeKey(i)
        self.Omega[target] = newSet
        self.updateBoundary(target)
        
    def add(self, S):
        """Add a (shallow) copy of a set S at an available index between
        0 and max + 1. Does not check for overlap! Calls updateBoundary,
        newIndex, and updates keyArchive."""
        if type(S) != set:
            raise ValueError, 'S must be set'
        index = self.newIndex()
        self.Omega[index] = S.copy()
        self.updateBoundary(index)
        self.keyArchive.add(index)
        
    def seed(self, seedCount = None, maxTry = 12):
        """seedCount defaults to 1"""
        """Add seedCount seed sets at available indices"""
        if seedCount == None:
            seedCount = 1
        else:
            if type(seedCount) != int or seedCount < 0:
                raise ValueError, "seedCount must be a positive integer"
        D = self.Domain.copy()
        (D.difference_update(self.completeOmega(i)) for i in self)
        numSeeds = 0
        for i in xrange(seedCount):
            for j in xrange(maxTry):
                p = utes.randomPoint(D)
                #print p
                S = utes.neighbors(p)
                #print S
                if S.issubset(D):
                    #print 'in if'
                    self.add(set([p]))
                    numSeeds += 1
                    #print 'prebreak, ', p, self
                    break
                    
        print "Created %s seeds out of %s requested" %(numSeeds, seedCount)
                            
    def indicator(self, index):
        M = np.zeros((self.n, self.n))
        for p in self[index]:
            M[p] = 1
        return M
        
    def bIndicator(self, index):
        M = np.zeros((self.n, self.n))
        for p in self[index,'b']:
            M[p] = 1
        return M
                        
    ##################################
    #                                #
    #  Graphing/Display Functions    #
    #                                #
    ##################################
    
    def display(self):
        return gUtes.returnPartition(self)
    
    def show(self):
        self.display()
    
    def __repr__(self):
        return "partition, n = %s, keys() = %s" % (self.n, self.keys())  
        
    #####################################
    #                                   #
    #   Matrix management functions     #
    #                                   #
    #####################################
    
    def globalOperator(self):
        return setup.globalOperator(self)
        
    def localOperator(self, index):
        indices =  [partition.indexing[a] for a in self[index]]
        return utes.csrSlice(partition.csrOp, indices)
                                
    def localOperator2(self, index):
        return utes.cooSlice3(partition.cooOp, [partition.indexing[a] for a in self[index]])
        
                      
            
if __name__ == '__main__':
    n = 100
    U = partition(n, 1./n)
    U.add(utes.graphball((20, 20), 10, U.Domain))
    U.add(utes.graphball((20,42), 10, U.Domain))
    U.add(utes.graphball((80,80), 10, U.Domain))
    W = partition(n, 1./n)
    W.seed(3)