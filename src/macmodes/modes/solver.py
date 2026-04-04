import os
import numpy as np
from dataclasses import dataclass
from petsc4py import PETSc
from slepc4py import SLEPc
from macmodes.modes.grids import get_grid, get_grid_params
from macmodes.modes.matrices import setup_mats_petsc
from macmodes.config import ModeConfig
from macmodes.constants import R_CMB, SECONDS_PER_YEAR, OMEGA

@dataclass
class Modes:
    nconv: int
    evals: np.ndarray[complex]
    evecs: np.ndarray[complex]
    freqs: np.ndarray
    damps: np.ndarray
    qs: np.ndarray
    pers_yrs: np.ndarray

def solve_slepc(A: PETSc.Mat, B: PETSc.Mat, n_eigs: int):
    ''' solves eigenvalue problem Ax +omega*Bx = 0'''

    print("solve_slepc: starting execution")
    vr, vi = A.createVecs() 

    E = SLEPc.EPS().create()
    E.setOperators(A, B)
    E.setType(SLEPc.EPS.Type.KRYLOVSCHUR) # set solver type
    E.setDimensions(n_eigs, PETSc.DECIDE) # set number of eigenvalues to compute
    E.setWhichEigenpairs(E.Which.SMALLEST_REAL) # look for smallest frequency
    E.setProblemType(SLEPc.EPS.ProblemType.GNHEP) # geneneralized non-hermitian

    E.solve()
    
    # print solver info
    
    print("")
    print("solve_slepc: solver information")
    its = E.getIterationNumber()
    print("Number of iterations of the method: %i" % its)
    sol_type = E.getType()
    print("Solution method: %s" % sol_type)
    nev, ncv, mpd = E.getDimensions()
    print("Number of requested eigenvalues: %i" % nev)
    tol, maxit = E.getTolerances()
    print("Stopping condition: tol=%.4g, maxit=%d" % (tol, maxit))
    nconv = E.getConverged()
    print("Number of converged eigenpairs: %d" % nconv)

    # unpack into numpy arrays

    if nconv > 0:
        vec_size = A.getSize()[0]
        print(f"solve_petsc: vec_size = {vec_size}")
        vals = np.empty(nconv, dtype=complex)
        vecs = np.empty([vec_size, nconv], dtype=complex)
        for i in range(nconv):
            eigval = E.getEigenpair(i, vr, vi)
            error = E.computeError(i)
            vals[i] = eigval
            vecs[:,i] = vr.getValues([range(vec_size)]) + 1j * vi.getValues([range(vec_size)])

    print("solve_slepc: exiting")
    return vals, vecs, nconv

def sort_modes(evals: np.ndarray[complex],
               evecs: np.ndarray[complex],
               nconv: int) -> Modes:
    print(nconv)
    # Sort via frequency (real part)
    inds = np.argsort(np.abs(evals.real))
    evals = evals[inds]
    evecs = evecs[:, inds]

    return  Modes(
            nconv = nconv,
            evals = evals,
            evecs = evecs,
            freqs = evals.real,
            damps = evals.imag,
            qs = 0.5*np.abs(evals.real/evals.imag),
            pers_yrs = 2*np.pi/(evals.real*SECONDS_PER_YEAR*OMEGA),
            )

    
def compute_modes(r: np.ndarray[float],
                  grid_params: tuple[np.ndarray[float]],
                  config: ModeConfig) -> Modes:
    '''computes mac modes with grid defined by r and grid params, and parameters 
    from config'''

    # create petsc matrices
    print("compute_modes: entering setup_mats_petsc")
    A, B = setup_mats_petsc(r, grid_params, config)

    # solve eigenvalue problem (TO DO: include delta iteration in solver routine)
    # perhaps set flag for if I want delta_iteration in config

    print("compute_modes: entering solve_slepc")
    evals, evecs, nconv = solve_slepc(A, B, config.n_modes)

    # sort and process the modes (lowest to highest frequency)

    print("compute_modes: entering sort_modes")
    modes: Modes
    modes = sort_modes(evals, evecs, nconv)       

    return modes
