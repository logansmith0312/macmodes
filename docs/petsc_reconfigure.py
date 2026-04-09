#!/usr/bin/python3
if __name__ == '__main__':
  import sys
  import os
  sys.path.insert(0, os.path.abspath('config'))
  import configure
  configure_options = [
    '--download-cmake',
    '--download-fblaslapack',
    '--download-metis',
    '--download-mpich',
    '--download-mumps',
    '--download-parmetis',
    '--download-scalapack',
    '--download-superlu_dist',
    '--with-scalar-type=complex',
    'PETSC_ARCH=arch-linux-c-debug',
  ]
  configure.petsc_configure(configure_options)
