############################################################
#
# EISCOR (EIgenSolvers based on unitary CORe transformations)
#
############################################################

LAPACK_LIBS := -llapack -lblas

# directory where the library will be installed
INSTALLDIR := $(HOME)

# fortran compiler and flags
FC := gfortran
FFLAGS := -O3 -std=f95

ifdef LAPACK_LIBS
  FFLAGS := $(FFLAGS) -DHAVE_LAPACK
endif

