from DFJMminimizer import *

Vname = 'V.dat'

l = 5
ticks = 20
L = l*ticks
n = L + 2

V_ = Load(Vname)
V = expandV(V_, l)

U = domain(n, 1./n)
U.Omega[0] = U.Domain.copy()
U.update()

M, indexing = Operator_to_cooMatrix(U, 0, V, U.meshsize)
Mdense = M.todense()