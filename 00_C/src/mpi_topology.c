#include "mpi_topology.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

int nprocs, myrank;
MPI_Comm mpi_world_cart;
int np_dim[NDIMS];
int period[NDIMS];
cart_comm_1d comm_1d_x, comm_1d_y, comm_1d_z;
boundary_comm comm_boundary;

void mpi_topology_create() {
    int remain[3];

    MPI_Cart_create(MPI_COMM_WORLD, NDIMS, np_dim, period, 0, &mpi_world_cart);

    // X-direction
    remain[0] = 1; remain[1] = 0; remain[2] = 0;
    MPI_Cart_sub(mpi_world_cart, remain, &comm_1d_x.mpi_comm);
    MPI_Comm_rank(comm_1d_x.mpi_comm, &comm_1d_x.myrank);
    MPI_Comm_size(comm_1d_x.mpi_comm, &comm_1d_x.nprocs);
    MPI_Cart_shift(comm_1d_x.mpi_comm, 0, 1, &comm_1d_x.west_rank, &comm_1d_x.east_rank);

    // Y-direction
    remain[0] = 0; remain[1] = 1; remain[2] = 0;
    MPI_Cart_sub(mpi_world_cart, remain, &comm_1d_y.mpi_comm);
    MPI_Comm_rank(comm_1d_y.mpi_comm, &comm_1d_y.myrank);
    MPI_Comm_size(comm_1d_y.mpi_comm, &comm_1d_y.nprocs);
    MPI_Cart_shift(comm_1d_y.mpi_comm, 0, 1, &comm_1d_y.west_rank, &comm_1d_y.east_rank);

    // Z-direction
    remain[0] = 0; remain[1] = 0; remain[2] = 1;
    MPI_Cart_sub(mpi_world_cart, remain, &comm_1d_z.mpi_comm);
    MPI_Comm_rank(comm_1d_z.mpi_comm, &comm_1d_z.myrank);
    MPI_Comm_size(comm_1d_z.mpi_comm, &comm_1d_z.nprocs);
    MPI_Cart_shift(comm_1d_z.mpi_comm, 0, 1, &comm_1d_z.west_rank, &comm_1d_z.east_rank);
}

void mpi_topology_destroy() {
    MPI_Comm_free(&mpi_world_cart);
}

void mpi_boundary_create() {
    int i, j, k;
    int coords[3];
    MPI_Group group_cart, group_boundary;
    int total = np_dim[0] * np_dim[1] * np_dim[2];
    int *rank = malloc(total * sizeof(int));
    int curr_rank, n_rank = 0;

    MPI_Comm_group(mpi_world_cart, &group_cart);

    for (i = 0; i < np_dim[0]; i++) {
        for (j = 0; j < np_dim[1]; j++) {
            for (k = 0; k < np_dim[2]; k++) {
                coords[0] = i;
                coords[1] = j;
                coords[2] = k;

                if ((!period[0] && (i == 0 || i == np_dim[0]-1)) ||
                    (!period[1] && (j == 0 || j == np_dim[1]-1)) ||
                    (!period[2] && (k == 0 || k == np_dim[2]-1))) {
                    MPI_Cart_rank(mpi_world_cart, coords, &curr_rank);
                    rank[n_rank++] = curr_rank;
                }
            }
        }
    }

    if (n_rank > 0) {
        MPI_Group_incl(group_cart, n_rank, rank, &group_boundary);
        MPI_Comm_create(mpi_world_cart, group_boundary, &comm_boundary.mpi_comm);
    } else {
        comm_boundary.mpi_comm = MPI_COMM_NULL;
    }

    comm_boundary.myrank = -1;
    comm_boundary.nprocs = -1;

    if (comm_boundary.mpi_comm != MPI_COMM_NULL) {
        MPI_Comm_rank(comm_boundary.mpi_comm, &comm_boundary.myrank);
        MPI_Comm_size(comm_boundary.mpi_comm, &comm_boundary.nprocs);
    }

    free(rank);
    MPI_Group_free(&group_cart);
    if (n_rank > 0) MPI_Group_free(&group_boundary);
}