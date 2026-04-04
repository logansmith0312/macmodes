import numpy as np
from scipy.sparse import csr_matrix
import petsc4py
from petsc4py import PETSc

x = np.ones(5)
A = np.diag(x) * 1j
A_sp = csr_matrix(A) 

A_petsc = PETSc.Mat().createAIJ(size=A_sp.size, csr=(A_sp.indptr, A_sp.indices, A_sp.data))
A_petsc.assemble()

A_test = A_petsc.getValues(range(5), range(5))
A_test2 = A_petsc.getArray()

print(type(A_test))
print(A_test)
print("test 2")
print(A_test2)
print(type(A_petsc))
