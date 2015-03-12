import numpy as np
import utilities as utes
from twoDPartition import partition
import harmonicSquare as setup
import DFJM
import graphingUtilities as gUtes
import time

iterations = setup.iterations

framerate = 30

n = setup.n
U = partition(n, 1./n)
S = U.Domain
of, tf, ot, tt, oh = int(n*.2), int(n*.4), int(1.*n/3), int(2.*n/3), int(.5*n)

C = set([(of, n-tf), (tf, n-of)])
C.update(set([(i, tt) for i in xrange(n-of, n+1)]))
for i in xrange(oh, n+1):
    C.update(set([(i, j) for j in xrange(ot, tf+1)]))
for i in xrange(tf, oh):
    C.update(set([(i, j) for j in xrange(of, tf+1)]))
    
S.difference_update(C)
U.Domain = S.copy()
U.updateBoundary(ofDomain = True)
g = DFJM.globalLandscape(U)

if __name__ == '__main__':
    t0 = time.time()
    timeDeltas = list()
    for i in xrange(1):
        U = partition(n, 1./n)
        seedCount = 60
        S.difference_update(C)
        U.Domain = S.copy()
        U.updateBoundary(ofDomain = True)
        U.seed(seedCount)
        data = DFJM.runSearch(U, makeMovie = True, iterations=iterations, ipf=5)
        movieName = 'FMnumsets'+str(seedCount) 
        movieName += 'Lin'+str(setup.Clinear)
        movieName += 'Sq'+str(setup.Csquare)
        movieName += 'Tinit'+str(setup.Tinit)
        movieName += 'trial'+str(i)
        print time.time()-t0 ,' seconds'
        
        gUtes.makeMovie(data['movie'], movieName,fps = 30)
        print 'movie made titled', movieName
        tlast = sum(timeDeltas)
        print time.time()-tlast ,' seconds'
        print time.time()-t0, 'seconds total'
        timeDeltas.append(time.time()-sum(timeDeltas))
        utes.picklesave(data, movieName+'.dat')
    