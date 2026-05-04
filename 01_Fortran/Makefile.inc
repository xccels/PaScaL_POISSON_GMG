###############################################################################
# Compile setting with intel compiler and intel MPI.
# Modify flags according to the system and environment
###############################################################################

#FC = mpiifort

#FLAG  = -O3 -fPIC -cpp
#FLAG += -g
#FLAG += -fp-model=precise 
#FLAG += -check bounds

###############################################################################
# Compile setting with nvfortran and non-intel MPI.
# Modify flags according to the system and environment
###############################################################################

FC = mpif90

FLAG  = -O3 -fPIC -cpp 
