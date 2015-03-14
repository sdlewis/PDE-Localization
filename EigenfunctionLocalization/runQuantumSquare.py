import numpy as np
import utilities as utes
from twoDPartition import partition
import setup
import DFJM
import graphingUtilities as gUtes
import time
import itertools as it


if __name__ == '__main__' and 0 == 1:
    n = setup.n
    U = partition(n, 1./n)
    ite = setup.iterations

    
    t0 = time.time()
    num_seeds = setup.seeds

    U.seed(num_seeds)
    K = U.keys().copy()
    for i in xrange(3):
        for k in K:
            if k in U:
                U.evolve(k, pgrow = 1)
                
    data = DFJM.runSearch(U, iterations = ite, makeMovie = True, printUpdate = 100, ipf = 2)
    movieName = 'qs,'+'numsets'+str(num_seeds) + 'Lin'+str(setup.Clinear)
    movieName += 'Sq'+str(setup.Csquare)
    movieName += 'iter'+str(setup.iterations)
    movieName += 'Tinit'+str(setup.Tinit)
    movieName += 'n'+str(n)
    movieName += 'trial0'
    t = time.time()-t0
    print '%s minutes, %s seconds' % (int(t)/60, t%60)
    #gUtes.makeMovie(data['movie'], movieName, fps = 10)
    utes.picklesave(data, movieName+'.dat')
    print 'data saved under ', movieName
    #print time.time()-t0 ,' seconds'