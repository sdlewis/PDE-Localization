import twoDPartition
import utilities as utes
import graphingUtilities as gUtes
import numpy as np
import itertools as it

def __format_eigen_data(dat,k):
    new_dat = list()
    for l in xrange(k):
        new_f = np.zeros((80,80))
        for i,j in it.product(xrange(100), repeat= 2):
            if i%5 != 4 and j%5 != 4:
                new_f[4*(i/5) + i%5, 4*(j/5) + j%5] = dat[l][i,j]
        new_dat.append(new_f)
    return new_dat 
    
display = lambda U, d: gUtes.showPartOverFxnList(U,d,4,3)

if __name__ == "__main__":
    partition_data = utes.pickleload("FMnumsets120Lin1e-06Sq0.0035Tinit1e-11trial0.dat")
    U = partition_data['best']
    eigen_data = utes.pickleload("fmeigendata.dat")
    eigen_data = __format_eigen_data(eigen_data[1], 11)
    
    