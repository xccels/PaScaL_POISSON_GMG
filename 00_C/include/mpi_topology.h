#ifndef MPI_TOPOLOGY_H
#define MPI_TOPOLOGY_H

#include <mpi.h>
#include <stdbool.h>
#include "global.h"

#define NDIMS 3

extern int myrank, nprocs;
extern MPI_Comm mpi_world_cart;
extern int np_dim[NDIMS];
extern int period[NDIMS];  

typedef struct {
    int myrank, nprocs;
    int west_rank, east_rank;
    MPI_Comm mpi_comm;
} cart_comm_1d;

extern cart_comm_1d comm_1d_x, comm_1d_y, comm_1d_z;

typedef struct {
    int myrank, nprocs;
    MPI_Comm mpi_comm;
} boundary_comm;

extern boundary_comm comm_boundary;

void mpi_topology_create(void);
void mpi_topology_destroy(void);
void mpi_boundary_create(void);

#endif // MPI_TOPOLOGY_H
