import os
import numpy as np
from time import perf_counter
from dataclasses import dataclass
from petsc4py import PETSc
from slepc4py import SLEPc
from macmodes.modes.grids import get_grid, get_grid_params
from macmodes.modes.matrices import setup_mats_petsc_p
from macmodes.config import ModeParams, ModeFlags
from macmodes.constants import SECONDS_PER_YEAR, OMEGA

@dataclass
class Modes:
    n_modes: int
    evals: np.ndarray[complex]
    evecs: np.ndarray[complex]
    freqs: np.ndarray
    damps: np.ndarray
    qs: np.ndarray
    pers_yrs: np.ndarray

def solve_slepc(A: PETSc.Mat, B: PETSc.Mat, n_eigs: int,
                st_target_per_yrs: float):
    ''' solves eigenvalue problem Ax +omega*Bx = 0'''

    rank = PETSc.COMM_WORLD.getRank()

    if rank == 0:
        print("solve_slepc: starting execution")

    st_target = 2*np.pi/(st_target_per_yrs*SECONDS_PER_YEAR*OMEGA) 
    vr, vi = A.createVecs() 

    F=SLEPc.ST().create()
    F.setType(SLEPc.ST.Type.SINVERT)
    F.setShift(st_target)
    
    ksp = F.getKSP()
    ksp.setType(PETSc.KSP.Type.PREONLY)

    pc = ksp.getPC()
    pc.setType(PETSc.PC.Type.LU)
    pc.setFactorSolverType("superlu_dist")
    pc.setFactorShift(PETSc.Mat.FactorShiftType.POSITIVE_DEFINITE, amount=1e-6)
    pc.setFactorPivot(1e-6)

    E = SLEPc.EPS().create()
    E.setST(F)
    E.setOperators(A, B)
    E.setType(SLEPc.EPS.Type.KRYLOVSCHUR) # set solver type
    ncv = min(4*n_eigs, A.getSize()[0] - 1) # krylov subspace size
    E.setDimensions(n_eigs, ncv) # set number of eigenvalues to compute
    E.setTolerances(tol=1e-7, max_it=1000)
    E.setWhichEigenpairs(E.Which.TARGET_REAL) 
    E.setTarget(st_target)
    E.setProblemType(SLEPc.EPS.ProblemType.GNHEP) # geneneralized non-hermitian
    E.setFromOptions()

    #print("||A|| =", A.norm())
    #print("||B|| =", B.norm())

    #print("A nnz:", A.getInfo()['nz_used'])
    #print("B nnz:", B.getInfo()['nz_used'])

    E.solve()
    
    # print solver info
    its = E.getIterationNumber()
    sol_type = E.getType()
    nev, _, _ = E.getDimensions()
    tol, maxit = E.getTolerances()
    nconv = E.getConverged()
    if rank == 0:    
        print("")
        print("solve_slepc: solver information")
        print("Number of iterations of the method: %i" % its)
        print("Solution method: %s" % sol_type)
        print("Number of requested eigenvalues: %i" % nev)
        print("Stopping condition: tol=%.4g, maxit=%d" % (tol, maxit))
        print("Number of converged eigenpairs: %d" % nconv)

    # unpack into numpy arrays
    vec_size = A.getSize()[0] 

    vals = None
    vecs = None
    n_modes = 0

    if nconv == 0:
        if rank == 0:
            print("No converged eigenpairs at this target.")
            print("")
            print("solve_slepc: exiting")
        return None, None, 0

    if nconv > 0:
        n_modes = min(n_eigs,nconv) 

        # set up scatterers and sequential vectors where the local parallel vectors
        # will be gathered to
        scatter_r, vr_seq = setup_scatter_to_root(vr)
        scatter_i, vi_seq = setup_scatter_to_root(vi)

        if rank == 0:
            vals = np.empty(n_modes, dtype=complex)
            vecs = np.empty([vec_size, n_modes], dtype=complex)

        for i in range(n_modes):
            eigval = E.getEigenpair(i, vr, vi)

            if rank == 0:
                vals[i] = eigval

            scatter_r.scatter(vr, vr_seq, addv=PETSc.InsertMode.INSERT_VALUES)
            scatter_i.scatter(vi, vi_seq, addv=PETSc.InsertMode.INSERT_VALUES)

            if rank == 0:
                vecs[:,i] = vr_seq.getArray() + 1j * vi_seq.getArray()

    if rank == 0:
        print("")
        print("solve_slepc: exiting")
    return vals, vecs, n_modes

def sort_modes(evals: np.ndarray[complex],
               evecs: np.ndarray[complex],
               n_modes: int) -> Modes:
    # Sort via frequency (real part)
    inds = np.argsort(np.abs(evals.real))
    evals = evals[inds]
    evecs = evecs[:, inds]

    return  Modes(
            n_modes = n_modes,
            evals = evals,
            evecs = evecs,
            freqs = evals.real,
            damps = evals.imag,
            qs = 0.5*np.abs(evals.real/evals.imag),
            pers_yrs = 2*np.pi/(evals.real*SECONDS_PER_YEAR*OMEGA),
            )

def setup_scatter_to_root(parallel_vec: PETSc.Vec):
    """
    Create a VecScatter that gathers a distributed PETSc Vec
    into a sequential Vec on rank 0.
    """
    comm = parallel_vec.getComm()
    rank = comm.getRank()
    global_size = parallel_vec.getSize()

    if rank == 0:
        # full size sequential vector on root 
        seq_vec = PETSc.Vec().createSeq(global_size, comm=PETSc.COMM_SELF)
    else:
        # size 0 vector on non-root ranks
        seq_vec = PETSc.Vec().createSeq(0, comm = PETSc.COMM_SELF)

    # PETSc scatter object 
    scatter, _, = PETSc.Scatter.toZero(parallel_vec)

    return scatter, seq_vec
    
def compute_modes(r: np.ndarray[float],
                  grid_params: tuple[np.ndarray[float]],
                  params: ModeParams,
                  flags: ModeFlags) -> Modes:
    '''computes mac modes with grid defined by r and grid params, and parameters 
    from config'''

    # create petsc matrices
    rank = PETSc.COMM_WORLD.getRank()
    if rank == 0:
        print("compute_modes: entering setup_mats_petsc")
    
    if rank == 0:
        tic = perf_counter()
    
    A, B = setup_mats_petsc_p(r, grid_params, params, flags)

    if rank == 0:
        tok = perf_counter()
        print(f"Matrix assembly time: {tok - tic:.3f} seconds")

    #print("compute_modes (DEBUG): matrix A")
    #print(f"size A: {A.getSize()}")
    #print(A.getValues(range(15),range(15)))
    #print("compute_modes (DEBUG): matrix B")
    #print(f"size B: {B.getSize()}")
    #print(B.getValues(range(15), range(15)))
    #print(f"A[10,13]: {A.getValue(10,13)}")

    #size = params.ell_max*params.n_rad*5//2
    #np.savez("matrices_slepc.npz", A=A.getValues(range(size),range(size,)), B=B.getValues(range(size),range(size)))

    # solve eigenvalue problem (TO DO: include delta iteration in solver routine)
    # perhaps set flag for if I want delta_iteration in config

    if rank == 0:
        print("compute_modes: entering solve_slepc")

    if rank == 0:
        tic = perf_counter()

    evals, evecs, n_modes = solve_slepc(A, B, params.n_modes,
                                        params.st_target_per_yrs)

    if rank == 0:
        tok = perf_counter()
        print(f"SLEPc solve time: {tok - tic:.3f} seconds")
    # sort and process the modes (lowest to highest frequency)

    if rank == 0:
        if n_modes > 0:
            print("compute_modes: entering sort_modes")
            modes: Modes
            modes = sort_modes(evals, evecs, n_modes)       
        else:
            print("compute_modes: no converged eigenpairs; returning empty modes")
            modes = Modes(
                n_modes = 0,
                evals = np.array([], dtype=complex),
                evecs = np.array((A.getSize()[0], 0), dtype=complex),
                freqs = np.array([]),
                damps = np.array([]),
                qs = np.array([]),
                pers_yrs = np.array([]),
                )

    else:
        modes = None

    return modes
