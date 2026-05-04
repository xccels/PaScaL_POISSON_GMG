# MPI C compiler
CC      = mpicc
CFLAGS  = -O3 -march=native -ffast-math -funroll-loops -fopenmp
LDFLAGS = -lm