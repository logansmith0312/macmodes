from petsc4py import PETSc

A = PETSc.Mat().createAIJ([4,4], comm=PETSc.COMM_WORLD)
A.setUp()

# simple diagonal matrix
Istart, Iend = A.getOwnershipRange()
for i in range(Istart, Iend):
    A.setValue(i, i, 2.0)
A.assemble()

ksp = PETSc.KSP().create()
ksp.setOperators(A)
ksp.setType('preonly')
ksp.getPC().setType('lu')

ksp.setFromOptions()
b = A.createVecRight()
x = A.createVecRight()
b.set(1.0)   # or any RHS
ksp.solve(b, x)

print("LU factorization succeeded")

