import numpy as np
from petsc4py import PETSc


A11 = PETSc.Mat().create()
n = 2
A11.setSizes([n,n])
A11.setType('aij')

for i in range(n):
    A11.setValue(i,i,1)
A11.assemble()

ell_max = 2
zero = np.zeros([n,n])
nest = [ [A11 , None], [None, A11]]

A = PETSc.Mat().createNest(nest)

A12 = PETSc.Mat().create()
n = 2
A12.setSizes([n,n])
A12.setType('aij')

for i in range(n):
    A12.setValue(i,i,4)
A12.assemble()

Alist = A.getNestSubMats()



A_np = A.getNestSubMatrix(1,1).getValues(range(n),range(n))
print(A_np)



