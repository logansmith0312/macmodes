from petsc4py import PETSc
print("Rank", PETSc.COMM_WORLD.getRank(), "of", PETSc.COMM_WORLD.getSize())
