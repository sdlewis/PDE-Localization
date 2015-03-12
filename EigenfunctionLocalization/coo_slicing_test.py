import itertools as it
from DFJM import *

if __name__=='__main__':
    n = setup.n
    U = partition(n, 1./n)
    U.add(utes.graphball((200,200),150, U.Domain))
    
    