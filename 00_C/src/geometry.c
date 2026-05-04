#include "geometry.h"
#include "mpi_topology.h"
#include "para_range.h"  // para_range assumed here
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>

void geometry_domain_create(domain *gdm, int nx, int ny, int nz,
                            double ox, double oy, double oz,
                            double lx, double ly, double lz,
                            double ax, double ay, double az,
                            int period[3])
{
    int i;

    gdm->nx = nx;
    gdm->ny = ny;
    gdm->nz = nz;
    gdm->ox = ox;
    gdm->oy = oy;
    gdm->oz = oz;
    gdm->lx = lx;
    gdm->ly = ly;
    gdm->lz = lz;
    for (i = 0; i < 3; i++) {
        gdm->is_periodic[i] = period[i];
    }

    gdm->dxm = (double *)calloc(gdm->nx+2, sizeof(double));
    gdm->dym = (double *)calloc(gdm->ny+2, sizeof(double));
    gdm->dzm = (double *)calloc(gdm->nz+2, sizeof(double));

    gdm->dxg = (double *)calloc(gdm->nx+2, sizeof(double));
    gdm->dyg = (double *)calloc(gdm->ny+2, sizeof(double));
    gdm->dzg = (double *)calloc(gdm->nz+2, sizeof(double));

    gdm->xg = (double *)calloc(gdm->nx+2, sizeof(double));
    gdm->yg = (double *)calloc(gdm->ny+2, sizeof(double));
    gdm->zg = (double *)calloc(gdm->nz+2, sizeof(double));

    if (ax == 1.0) 
    {
        gdm->dxm[0] = 0.0;
        gdm->dxg[0] = 0.0;
        gdm->xg[0] = gdm->ox;
        for (i = 1; i <= gdm->nx; i++) {
            gdm->dxm[i] = gdm->lx / gdm->nx;
            gdm->dxg[i] = 0.5 * (gdm->dxm[i] + gdm->dxm[i-1]);
            gdm->xg[i] = gdm->xg[i-1] + gdm->dxg[i];
        }
        gdm->dxm[gdm->nx+1] = 0.0;
        gdm->dxg[gdm->nx+1] = 0.5 * (gdm->dxm[gdm->nx+1] + gdm->dxm[gdm->nx]);
        gdm->xg[gdm->nx+1] = gdm->xg[gdm->nx] + gdm->dxg[gdm->nx+1];
    } 
    else 
    {
        gdm->dxm[0] = 0.0;
        gdm->dxm[1] = 0.5 * lx * (ax - 1.0) / (pow(ax, gdm->nx/2) - 1.0);
        for (i = 2; i <= gdm->nx/2; i++) 
        {
            gdm->dxm[i] = gdm->dxm[i-1] * ax;
            gdm->dxm[gdm->nx - i + 1] = gdm->dxm[i];
        }
        gdm->dxm[gdm->nx] = gdm->dxm[1];
        gdm->dxm[gdm->nx+1] = 0.0;

        gdm->dxg[0] = 0.0;
        gdm->xg[0] = gdm->ox;
        for (i = 1; i <= gdm->nx; i++) {
            gdm->dxg[i] = 0.5 * (gdm->dxm[i] + gdm->dxm[i-1]);
            gdm->xg[i] = gdm->xg[i-1] + gdm->dxg[i];
        }
        gdm->dxg[gdm->nx+1] = 0.5 * (gdm->dxm[gdm->nx+1] + gdm->dxm[gdm->nx]);
        gdm->xg[gdm->nx+1] = gdm->xg[gdm->nx] + gdm->dxg[gdm->nx+1];
    }

    if (ay == 1.0) 
    {
        gdm->dym[0] = 0.0;
        gdm->dyg[0] = 0.0;
        gdm->yg[0] = gdm->oy;
        for (i = 1; i <= gdm->ny; i++) {
            gdm->dym[i] = gdm->ly / gdm->ny;
            gdm->dyg[i] = 0.5 * (gdm->dym[i] + gdm->dym[i-1]);
            gdm->yg[i] = gdm->yg[i-1] + gdm->dyg[i];
        }
        gdm->dym[gdm->ny+1] = 0.0;
        gdm->dyg[gdm->ny+1] = 0.5 * (gdm->dym[gdm->ny+1] + gdm->dym[gdm->ny]);
        gdm->yg[gdm->ny+1] = gdm->yg[gdm->ny] + gdm->dyg[gdm->ny+1];
    } 
    else 
    {
        gdm->dym[0] = 0.0;
        gdm->dym[1] = 0.5 * ly * (ay - 1.0) / (pow(ay, gdm->ny/2) - 1.0);
        for (i = 2; i <= gdm->ny/2; i++) 
        {
            gdm->dym[i] = gdm->dym[i-1] * ay;
            gdm->dym[gdm->ny - i + 1] = gdm->dym[i];
        }
        gdm->dym[gdm->ny] = gdm->dym[1];
        gdm->dym[gdm->ny+1] = 0.0;

        gdm->dyg[0] = 0.0;
        gdm->yg[0] = gdm->oy;
        for (i = 1; i <= gdm->ny; i++) {
            gdm->dyg[i] = 0.5 * (gdm->dym[i] + gdm->dym[i-1]);
            gdm->yg[i] = gdm->yg[i-1] + gdm->dyg[i];
        }
        gdm->dyg[gdm->ny+1] = 0.5 * (gdm->dym[gdm->ny+1] + gdm->dym[gdm->ny]);
        gdm->yg[gdm->ny+1] = gdm->yg[gdm->ny] + gdm->dyg[gdm->ny+1];
    }

    if (az == 1.0) 
    {
        gdm->dzm[0] = 0.0;
        gdm->dzg[0] = 0.0;
        gdm->zg[0] = gdm->oz;
        for (i = 1; i <= gdm->nz; i++) {
            gdm->dzm[i] = gdm->lz / gdm->nz;
            gdm->dzg[i] = 0.5 * (gdm->dzm[i] + gdm->dzm[i-1]);
            gdm->zg[i] = gdm->zg[i-1] + gdm->dzg[i];
        }
        gdm->dzm[gdm->nz+1] = 0.0;
        gdm->dzg[gdm->nz+1] = 0.5 * (gdm->dzm[gdm->nz+1] + gdm->dzm[gdm->nz]);
        gdm->zg[gdm->nz+1] = gdm->zg[gdm->nz] + gdm->dzg[gdm->nz+1];
    } 
    else 
    {
        gdm->dzm[0] = 0.0;
        gdm->dzm[1] = 0.5 * lz * (az - 1.0) / (pow(az, gdm->nz/2) - 1.0);
        for (i = 2; i <= gdm->nz/2; i++) 
        {
            gdm->dzm[i] = gdm->dzm[i-1] * az;
            gdm->dzm[gdm->nz - i + 1] = gdm->dzm[i];
        }
        gdm->dzm[gdm->nz] = gdm->dzm[1];
        gdm->dzm[gdm->nz+1] = 0.0;

        gdm->dzg[0] = 0.0;
        gdm->zg[0] = gdm->oz;
        for (i = 1; i <= gdm->nz; i++) {
            gdm->dzg[i] = 0.5 * (gdm->dzm[i] + gdm->dzm[i-1]);
            gdm->zg[i] = gdm->zg[i-1] + gdm->dzg[i];
        }
        gdm->dzg[gdm->nz+1] = 0.5 * (gdm->dzm[gdm->nz+1] + gdm->dzm[gdm->nz]);
        gdm->zg[gdm->nz+1] = gdm->zg[gdm->nz] + gdm->dzg[gdm->nz+1];
    }

}


void geometry_domain_destroy(domain *gdm) 
{
    free(gdm->dxm);
    free(gdm->dym);
    free(gdm->dzm);
    free(gdm->dxg);
    free(gdm->dyg);
    free(gdm->dzg);
    free(gdm->xg);
    free(gdm->yg);
    free(gdm->zg);

    gdm->dxm = gdm->dym = gdm->dzm = NULL;
    gdm->dxg = gdm->dyg = gdm->dzg = NULL;
    gdm->xg = gdm->yg = gdm->zg = NULL;
}

void geometry_subdomain_create(subdomain *sdm, const domain *gdm) {
    // extern cart_comm_1d comm_1d_x, comm_1d_y, comm_1d_z;

    for (int i = 0; i < 3; i++) {
       sdm->is_periodic[i] = gdm->is_periodic[i];
    }

    if(comm_1d_x.nprocs == 1) {
        sdm->is_aggregated[0] = 1;  
    } else {
        sdm->is_aggregated[0] = 0; 
    }
    if(comm_1d_y.nprocs == 1) {
        sdm->is_aggregated[1] = 1;  
    } else {
        sdm->is_aggregated[1] = 0; 
    }
    if(comm_1d_z.nprocs == 1) {
        sdm->is_aggregated[2] = 1;  
    } else {
        sdm->is_aggregated[2] = 0; 
    }

    para_range(1, gdm->nx, comm_1d_x.nprocs, comm_1d_x.myrank, &sdm->ista, &sdm->iend);
    sdm->nx = sdm->iend - sdm->ista + 1;
    para_range(1, gdm->ny, comm_1d_y.nprocs, comm_1d_y.myrank, &sdm->jsta, &sdm->jend);
    sdm->ny = sdm->jend - sdm->jsta + 1;
    para_range(1, gdm->nz, comm_1d_z.nprocs, comm_1d_z.myrank, &sdm->ksta, &sdm->kend);
    sdm->nz = sdm->kend - sdm->ksta + 1;

    

    int npx = sdm->nx + 2, npy = sdm->ny + 2, npz = sdm->nz + 2;

    sdm->dxm = malloc(npx * sizeof(double));
    sdm->dym = malloc(npy * sizeof(double));
    sdm->dzm = malloc(npz * sizeof(double));
    sdm->dxg = malloc(npx * sizeof(double));
    sdm->dyg = malloc(npy * sizeof(double));
    sdm->dzg = malloc(npz * sizeof(double));
    sdm->xg  = malloc(npx * sizeof(double));
    sdm->yg  = malloc(npy * sizeof(double));
    sdm->zg  = malloc(npz * sizeof(double));

    for (int i = 0; i < npx; i++)
    {
        sdm->dxm[i] = gdm->dxm[i + sdm->ista - 1];
        sdm->dxg[i] = gdm->dxg[i + sdm->ista - 1];
        sdm->xg[i]  = gdm->xg[i + sdm->ista - 1];
    }
    for (int j = 0; j < npy; j++)
    {
        sdm->dym[j] = gdm->dym[j + sdm->jsta - 1];
        sdm->dyg[j] = gdm->dyg[j + sdm->jsta - 1];
        sdm->yg[j]  = gdm->yg[j + sdm->jsta - 1];
    }
    for (int k = 0; k < npz; k++)
    {
        sdm->dzm[k] = gdm->dzm[k + sdm->ksta - 1];
        sdm->dzg[k] = gdm->dzg[k + sdm->ksta - 1];
        sdm->zg[k]  = gdm->zg[k + sdm->ksta - 1];
    }
    

    // memcpy(sdm->dxm, gdm->dxm + sdm->ista - 1, npx * sizeof(double));
    // memcpy(sdm->dym, gdm->dym + sdm->jsta - 1, npy * sizeof(double));
    // memcpy(sdm->dzm, gdm->dzm + sdm->ksta - 1, npz * sizeof(double));
    // memcpy(sdm->dxg, gdm->dxg + sdm->ista - 1, npx * sizeof(double));
    // memcpy(sdm->dyg, gdm->dyg + sdm->jsta - 1, npy * sizeof(double));
    // memcpy(sdm->dzg, gdm->dzg + sdm->ksta - 1, npz * sizeof(double));
    // memcpy(sdm->xg,  gdm->xg + sdm->ista - 1, npx * sizeof(double));
    // memcpy(sdm->yg,  gdm->yg + sdm->jsta - 1, npy * sizeof(double));
    // memcpy(sdm->zg,  gdm->zg + sdm->ksta - 1, npz * sizeof(double));

    sdm->ox = sdm->xg[1] - 0.5 * sdm->dxm[1];
    sdm->oy = sdm->yg[1] - 0.5 * sdm->dym[1];
    sdm->oz = sdm->zg[1] - 0.5 * sdm->dzm[1];

    size_t size3d = npx * npy * npz;
    sdm->x = calloc(size3d, sizeof(double));
    sdm->b = calloc(size3d, sizeof(double));
    sdm->r = calloc(size3d, sizeof(double));

    sdm->is_x0_boundary = (comm_1d_x.west_rank == MPI_PROC_NULL);
    sdm->is_x1_boundary = (comm_1d_x.east_rank == MPI_PROC_NULL);
    sdm->is_y0_boundary = (comm_1d_y.west_rank == MPI_PROC_NULL);
    sdm->is_y1_boundary = (comm_1d_y.east_rank == MPI_PROC_NULL);
    sdm->is_z0_boundary = (comm_1d_z.west_rank == MPI_PROC_NULL);
    sdm->is_z1_boundary = (comm_1d_z.east_rank == MPI_PROC_NULL);
}

void geometry_subdomain_destroy(subdomain *sdm) {
    free(sdm->dxm); free(sdm->dym); free(sdm->dzm);
    free(sdm->dxg); free(sdm->dyg); free(sdm->dzg);
    free(sdm->xg);  free(sdm->yg);  free(sdm->zg);
    free(sdm->x);   free(sdm->b);   free(sdm->r);
}

void geometry_subdomain_ddt_create(subdomain *sdm)
{
    int sizes[3], subsizes[3], starts[3];

    sizes[0] = sdm->nx + 2;
    sizes[1] = sdm->ny + 2;
    sizes[2] = sdm->nz + 2;

    // Inner domain
    subsizes[0] = sdm->nx;
    subsizes[1] = sdm->ny;
    subsizes[2] = sdm->nz;
    starts[0] = 1;
    starts[1] = 1;
    starts[2] = 1;
    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &sdm->ddt_inner_domain);
    MPI_Type_commit(&sdm->ddt_inner_domain);

    // yz_plane_xn
    subsizes[0] = 1; subsizes[1] = sdm->ny + 2; subsizes[2] = sdm->nz + 2;
    starts[0] = sdm->nx; starts[1] = 0; starts[2] = 0;
    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &sdm->ddt_yz_plane_xn);
    MPI_Type_commit(&sdm->ddt_yz_plane_xn);

    // yz_plane_x0
    starts[0] = 0;
    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &sdm->ddt_yz_plane_x0);
    MPI_Type_commit(&sdm->ddt_yz_plane_x0);

    // yz_plane_x1
    starts[0] = 1;
    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &sdm->ddt_yz_plane_x1);
    MPI_Type_commit(&sdm->ddt_yz_plane_x1);

    // yz_plane_xn1
    starts[0] = sdm->nx + 1;
    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &sdm->ddt_yz_plane_xn1);
    MPI_Type_commit(&sdm->ddt_yz_plane_xn1);

    // xz_plane_yn
    subsizes[0] = sdm->nx + 2; subsizes[1] = 1; subsizes[2] = sdm->nz + 2;
    starts[0] = 0; starts[1] = sdm->ny; starts[2] = 0;
    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &sdm->ddt_xz_plane_yn);
    MPI_Type_commit(&sdm->ddt_xz_plane_yn);

    // xz_plane_y0
    starts[1] = 0;
    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &sdm->ddt_xz_plane_y0);
    MPI_Type_commit(&sdm->ddt_xz_plane_y0);

    // xz_plane_y1
    starts[1] = 1;
    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &sdm->ddt_xz_plane_y1);
    MPI_Type_commit(&sdm->ddt_xz_plane_y1);

    // xz_plane_yn1
    starts[1] = sdm->ny + 1;
    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &sdm->ddt_xz_plane_yn1);
    MPI_Type_commit(&sdm->ddt_xz_plane_yn1);

    // xy_plane_zn
    subsizes[1] = sdm->ny + 2; subsizes[2] = 1;
    starts[1] = 0; starts[2] = sdm->nz;
    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &sdm->ddt_xy_plane_zn);
    MPI_Type_commit(&sdm->ddt_xy_plane_zn);

    // xy_plane_z0
    starts[2] = 0;
    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &sdm->ddt_xy_plane_z0);
    MPI_Type_commit(&sdm->ddt_xy_plane_z0);

    // xy_plane_z1
    starts[2] = 1;
    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &sdm->ddt_xy_plane_z1);
    MPI_Type_commit(&sdm->ddt_xy_plane_z1);

    // xy_plane_zn1
    starts[2] = sdm->nz + 1;
    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &sdm->ddt_xy_plane_zn1);
    MPI_Type_commit(&sdm->ddt_xy_plane_zn1);
}

void geometry_subdomain_ddt_destroy(subdomain *sdm) {
    int ierr;
    MPI_Type_free(&sdm->ddt_yz_plane_x0);
    MPI_Type_free(&sdm->ddt_yz_plane_x1);
    MPI_Type_free(&sdm->ddt_yz_plane_xn);
    MPI_Type_free(&sdm->ddt_yz_plane_xn1);

    MPI_Type_free(&sdm->ddt_xz_plane_y0);
    MPI_Type_free(&sdm->ddt_xz_plane_y1);
    MPI_Type_free(&sdm->ddt_xz_plane_yn);
    MPI_Type_free(&sdm->ddt_xz_plane_yn1);

    MPI_Type_free(&sdm->ddt_xy_plane_z0);
    MPI_Type_free(&sdm->ddt_xy_plane_z1);
    MPI_Type_free(&sdm->ddt_xy_plane_zn);
    MPI_Type_free(&sdm->ddt_xy_plane_zn1);
}

// void geometry_halocell_update(double *u, subdomain *sdm) {
//     extern cart_comm_1d comm_1d_x, comm_1d_y, comm_1d_z;

//     MPI_Request request[4];

//     // double *U = &u[0][0][0];

//     // // x-direction
//     // MPI_Isend(U, 1, sdm->ddt_yz_plane_xn,  comm_1d_x.east_rank, 111, comm_1d_x.mpi_comm, &request[0]);
//     // MPI_Irecv(U, 1, sdm->ddt_yz_plane_x0,  comm_1d_x.west_rank, 111, comm_1d_x.mpi_comm, &request[1]);
//     // MPI_Isend(U, 1, sdm->ddt_yz_plane_x1,  comm_1d_x.west_rank, 222, comm_1d_x.mpi_comm, &request[2]);
//     // MPI_Irecv(U, 1, sdm->ddt_yz_plane_xn1, comm_1d_x.east_rank, 222, comm_1d_x.mpi_comm, &request[3]);
//     // MPI_Waitall(4, request, MPI_STATUSES_IGNORE);

//     // // y-direction
//     // MPI_Isend(U, 1, sdm->ddt_xz_plane_yn,  comm_1d_y.east_rank, 333, comm_1d_y.mpi_comm, &request[0]);
//     // MPI_Irecv(U, 1, sdm->ddt_xz_plane_y0,  comm_1d_y.west_rank, 333, comm_1d_y.mpi_comm, &request[1]);
//     // MPI_Isend(U, 1, sdm->ddt_xz_plane_y1,  comm_1d_y.west_rank, 444, comm_1d_y.mpi_comm, &request[2]);
//     // MPI_Irecv(U, 1, sdm->ddt_xz_plane_yn1, comm_1d_y.east_rank, 444, comm_1d_y.mpi_comm, &request[3]);
//     // MPI_Waitall(4, request, MPI_STATUSES_IGNORE);

//     // // z-direction
//     // MPI_Isend(U, 1, sdm->ddt_xy_plane_zn,  comm_1d_z.east_rank, 555, comm_1d_z.mpi_comm, &request[0]);
//     // MPI_Irecv(U, 1, sdm->ddt_xy_plane_z0,  comm_1d_z.west_rank, 555, comm_1d_z.mpi_comm, &request[1]);
//     // MPI_Isend(U, 1, sdm->ddt_xy_plane_z1,  comm_1d_z.west_rank, 666, comm_1d_z.mpi_comm, &request[2]);
//     // MPI_Irecv(U, 1, sdm->ddt_xy_plane_zn1, comm_1d_z.east_rank, 666, comm_1d_z.mpi_comm, &request[3]);
//     // MPI_Waitall(4, request, MPI_STATUSES_IGNORE);

//     // x-direction
//     MPI_Isend(u, 1, sdm->ddt_yz_plane_xn,  comm_1d_x.east_rank, 111, comm_1d_x.mpi_comm, &request[0]);
//     MPI_Irecv(u, 1, sdm->ddt_yz_plane_x0,  comm_1d_x.west_rank, 111, comm_1d_x.mpi_comm, &request[1]);
//     MPI_Isend(u, 1, sdm->ddt_yz_plane_x1,  comm_1d_x.west_rank, 222, comm_1d_x.mpi_comm, &request[2]);
//     MPI_Irecv(u, 1, sdm->ddt_yz_plane_xn1, comm_1d_x.east_rank, 222, comm_1d_x.mpi_comm, &request[3]);
//     MPI_Waitall(4, request, MPI_STATUSES_IGNORE);

//     // y-direction
//     MPI_Isend(u, 1, sdm->ddt_xz_plane_yn,  comm_1d_y.east_rank, 333, comm_1d_y.mpi_comm, &request[0]);
//     MPI_Irecv(u, 1, sdm->ddt_xz_plane_y0,  comm_1d_y.west_rank, 333, comm_1d_y.mpi_comm, &request[1]);
//     MPI_Isend(u, 1, sdm->ddt_xz_plane_y1,  comm_1d_y.west_rank, 444, comm_1d_y.mpi_comm, &request[2]);
//     MPI_Irecv(u, 1, sdm->ddt_xz_plane_yn1, comm_1d_y.east_rank, 444, comm_1d_y.mpi_comm, &request[3]);
//     MPI_Waitall(4, request, MPI_STATUSES_IGNORE);

//     // z-direction
//     MPI_Isend(u, 1, sdm->ddt_xy_plane_zn,  comm_1d_z.east_rank, 555, comm_1d_z.mpi_comm, &request[0]);
//     MPI_Irecv(u, 1, sdm->ddt_xy_plane_z0,  comm_1d_z.west_rank, 555, comm_1d_z.mpi_comm, &request[1]);
//     MPI_Isend(u, 1, sdm->ddt_xy_plane_z1,  comm_1d_z.west_rank, 666, comm_1d_z.mpi_comm, &request[2]);
//     MPI_Irecv(u, 1, sdm->ddt_xy_plane_zn1, comm_1d_z.east_rank, 666, comm_1d_z.mpi_comm, &request[3]);
//     MPI_Waitall(4, request, MPI_STATUSES_IGNORE);
// }

void geometry_halocell_update_selectively(double *u, subdomain *sdm, int *is_serial) {

    MPI_Request request[4];

    if (!is_serial[0]) {
        MPI_Isend(u, 1, sdm->ddt_yz_plane_xn,  comm_1d_x.east_rank, 111, comm_1d_x.mpi_comm, &request[0]);
        MPI_Irecv(u, 1, sdm->ddt_yz_plane_x0,  comm_1d_x.west_rank, 111, comm_1d_x.mpi_comm, &request[1]);
        MPI_Isend(u, 1, sdm->ddt_yz_plane_x1,  comm_1d_x.west_rank, 222, comm_1d_x.mpi_comm, &request[2]);
        MPI_Irecv(u, 1, sdm->ddt_yz_plane_xn1, comm_1d_x.east_rank, 222, comm_1d_x.mpi_comm, &request[3]);
        MPI_Waitall(4, request, MPI_STATUSES_IGNORE);
    }

    if (!is_serial[1]) {
        MPI_Isend(u, 1, sdm->ddt_xz_plane_yn,  comm_1d_y.east_rank, 333, comm_1d_y.mpi_comm, &request[0]);
        MPI_Irecv(u, 1, sdm->ddt_xz_plane_y0,  comm_1d_y.west_rank, 333, comm_1d_y.mpi_comm, &request[1]);
        MPI_Isend(u, 1, sdm->ddt_xz_plane_y1,  comm_1d_y.west_rank, 444, comm_1d_y.mpi_comm, &request[2]);
        MPI_Irecv(u, 1, sdm->ddt_xz_plane_yn1, comm_1d_y.east_rank, 444, comm_1d_y.mpi_comm, &request[3]);
        MPI_Waitall(4, request, MPI_STATUSES_IGNORE);
    }

    if (!is_serial[2]) {
        MPI_Isend(u, 1, sdm->ddt_xy_plane_zn,  comm_1d_z.east_rank, 555, comm_1d_z.mpi_comm, &request[0]);
        MPI_Irecv(u, 1, sdm->ddt_xy_plane_z0,  comm_1d_z.west_rank, 555, comm_1d_z.mpi_comm, &request[1]);
        MPI_Isend(u, 1, sdm->ddt_xy_plane_z1,  comm_1d_z.west_rank, 666, comm_1d_z.mpi_comm, &request[2]);
        MPI_Irecv(u, 1, sdm->ddt_xy_plane_zn1, comm_1d_z.east_rank, 666, comm_1d_z.mpi_comm, &request[3]);
        MPI_Waitall(4, request, MPI_STATUSES_IGNORE);
    }
}

// void geometry_boundary_value_calculate(double *u, subdomain *sdm, domain *gdm, boundary_comm *comm_boundary) {
//     const double PI = 3.14159265358979323846;
//     int nxs = sdm->nx, nys = sdm->ny, nzs = sdm->nz;
//     int nxg = gdm->nx, nyg = gdm->ny, nzg = gdm->nz;
//     double dudx0 = 0.0, dudx1 = 0.0, dudy0 = 0.0, dudy1 = 0.0, dudz0 = 0.0, dudz1 = 0.0;

//     // double ***bval_x = (double ***)malloc((nyg) * sizeof(double **));
//     // double ***bval_y = (double ***)malloc((nxg) * sizeof(double **));
//     // double ***bval_z = (double ***)malloc((nxg) * sizeof(double **));
//     // for (int j = 0; j < nyg; ++j) {
//     //     bval_x[j] = (double **)malloc((nzg) * sizeof(double *));
//     //     for (int k = 0; k < nzg; ++k) {
//     //         bval_x[j][k] = (double *)calloc(2, sizeof(double)); // 0和1边界
//     //     }
//     // }
//     // for (int i = 0; i < nxg; ++i) {
//     //     bval_y[i] = (double **)malloc((nzg) * sizeof(double *));
//     //     bval_z[i] = (double **)malloc((nyg) * sizeof(double *));
//     //     for (int k = 0; k < nzg; ++k) {
//     //         bval_y[i][k] = (double *)calloc(2, sizeof(double));
//     //     }
//     //     for (int j = 0; j < nyg; ++j) {
//     //         bval_z[i][j] = (double *)calloc(2, sizeof(double));
//     //     }
//     // }
//     double *bval_x = calloc((nyg) * (nzg) * 2, sizeof(double));
//     double *bval_y = calloc((nxg) * (nzg) * 2, sizeof(double));
//     double *bval_z = calloc((nxg) * (nyg) * 2, sizeof(double));


//     #ifndef INPLACE
//     // double ***bval_xg = (double ***)malloc((nyg) * sizeof(double **));
//     // double ***bval_yg = (double ***)malloc((nxg) * sizeof(double **));
//     // double ***bval_zg = (double ***)malloc((nxg) * sizeof(double **));
//     // for (int j = 0; j < nyg; ++j) {
//     //     bval_xg[j] = (double **)malloc((nzg) * sizeof(double *));
//     //     for (int k = 0; k < nzg; ++k) {
//     //         bval_xg[j][k] = (double *)calloc(2, sizeof(double)); // 0和1边界
//     //     }
//     // }
//     // for (int i = 0; i < nxg; ++i) {
//     //     bval_yg[i] = (double **)malloc((nzg) * sizeof(double *));
//     //     bval_zg[i] = (double **)malloc((nyg) * sizeof(double *));
//     //     for (int k = 0; k < nzg; ++k) {
//     //         bval_yg[i][k] = (double *)calloc(2, sizeof(double));
//     //     }
//     //     for (int j = 0; j < nyg; ++j) {
//     //         bval_zg[i][j] = (double *)calloc(2, sizeof(double));
//     //     }
//     // }
//     double *bval_xg = calloc((nyg) * (nzg) * 2, sizeof(double));
//     double *bval_yg = calloc((nxg) * (nzg) * 2, sizeof(double));
//     double *bval_zg = calloc((nxg) * (nyg) * 2, sizeof(double));

//     #endif

//     double gf_coeff = 0.25 / PI;

//     // Contributions from boundaries in x-direction
//     if (sdm->is_x0_boundary) {
//         double dxg1 = sdm->dxg[1];
//         double dxg2 = sdm->dxg[2];
//         double alpha = dxg2 / dxg1;
//         double w1 = gf_coeff / dxg1 * (2.0 + alpha) / (1.0 + alpha);
//         double w2 = gf_coeff / dxg2 / (1.0 + alpha);

//         for (int kp = 0; kp < nzs; ++kp) {
//             for (int jp = 0; jp < nys; ++jp) {
//                 int jrp = jp + sdm->jsta - 1;
//                 int krp = kp + sdm->ksta - 1;
                
//                 double dudx0 = -((u[1][jp][kp] - u[0][jp][kp]) * w1 - (u[2][jp][kp] - u[1][jp][kp]) * w2) * sdm->dym[jp] * sdm->dzm[kp];

//                 // Boundary treatment in x-direction
//                 for (int k = 0; k < nzg; ++k) {
//                     for (int j = 0; j < nyg; ++j) {
//                         // Lower boundary in x-direction x=0
//                         double drx = gdm->xg[0] - gdm->xg[1];
//                         double dry = gdm->yg[j] - gdm->yg[jrp];
//                         double drz = gdm->zg[k] - gdm->zg[krp];
//                         double dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_x[j][k][0] += dudx0 / dr;

//                         // Upper boundary in x-direction x=nx+1
//                         drx = gdm->xg[nxg + 1] - gdm->xg[1];
//                         dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_x[j][k][1] += dudx0 / dr;
//                     }
//                 }

//                 // Boundary treatment in y-direction
//                 for (int k = 0; k < nzg; ++k) {
//                     for (int i = 0; i < nxg; ++i) {
//                         // Lower boundary in y-direction
//                         double drx = gdm->xg[i] - gdm->xg[1];
//                         double dry = gdm->yg[0] - gdm->yg[jrp];
//                         double drz = gdm->zg[k] - gdm->zg[krp];
//                         double dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_y[i][k][0] += dudx0 / dr;

//                         // Upper boundary in y-direction
//                         dry = gdm->yg[nyg + 1] - gdm->yg[jrp];
//                         dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_y[i][k][1] += dudx0 / dr;
//                     }
//                 }

//                 // Boundary treatment in z-direction
//                 for (int j = 0; j < nyg; ++j) {
//                     for (int i = 0; i < nxg; ++i) {
//                         // Lower boundary in z-direction
//                         double drx = gdm->xg[i] - gdm->xg[1];
//                         double dry = gdm->yg[j] - gdm->yg[jrp];
//                         double drz = gdm->zg[0] - gdm->zg[krp];
//                         double dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_z[i][j][0] += dudx0 / dr;

//                         // Upper boundary in z-direction
//                         drz = gdm->zg[nzg + 1] - gdm->zg[krp];
//                         dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_z[i][j][1] += dudx0 / dr;
//                     }
//                 }
//             }
//         }
//     }

//     // Contributions from boundaries in x-direction
//     if (sdm->is_x1_boundary) {
//         double dxg1 = sdm->dxg[nxs-2];
//         double dxg2 = sdm->dxg[nxs-1];
//         double alpha = dxg1 / dxg2;
//         double w1 = gf_coeff / dxg2 * (2.0 + alpha) / (1.0 + alpha);
//         double w2 = gf_coeff / dxg1 / (1.0 + alpha);

//         for (int kp = 0; kp < nzs; ++kp) {
//             for (int jp = 0; jp < nys; ++jp) {
//                 int jrp = jp + sdm->jsta - 1;
//                 int krp = kp + sdm->ksta - 1;
                
//                 double dudx1 = -((u[nxs-1][jp][kp] - u[nxs-2][jp][kp]) * w1 - (u[nxs-2][jp][kp] - u[nxs-3][jp][kp]) * w2) * sdm->dym[jp] * sdm->dzm[kp];

//                 // Boundary treatment in x-direction
//                 for (int k = 0; k < nzg; ++k) {
//                     for (int j = 0; j < nyg; ++j) {
//                         // Lower boundary in x-direction x=0
//                         double drx = gdm->xg[0] - gdm->xg[nxg];
//                         double dry = gdm->yg[j] - gdm->yg[jrp];
//                         double drz = gdm->zg[k] - gdm->zg[krp];
//                         double dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_x[j][k][0] += dudx1 / dr;

//                         // Upper boundary in x-direction x=nx+1
//                         drx = gdm->xg[nxg + 1] - gdm->xg[nxg];
//                         dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_x[j][k][1] += dudx1 / dr;
//                     }
//                 }

//                 // Boundary treatment in y-direction
//                 for (int k = 0; k < nzg; ++k) {
//                     for (int i = 0; i < nxg; ++i) {
//                         // Lower boundary in y-direction
//                         double drx = gdm->xg[i] - gdm->xg[nxg];
//                         double dry = gdm->yg[0] - gdm->yg[jrp];
//                         double drz = gdm->zg[k] - gdm->zg[krp];
//                         double dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_y[i][k][0] += dudx1 / dr;

//                         // Upper boundary in y-direction
//                         dry = gdm->yg[nyg + 1] - gdm->yg[jrp];
//                         dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_y[i][k][1] += dudx1 / dr;
//                     }
//                 }

//                 // Boundary treatment in z-direction
//                 for (int j = 0; j < nyg; ++j) {
//                     for (int i = 0; i < nxg; ++i) {
//                         // Lower boundary in z-direction
//                         double drx = gdm->xg[i] - gdm->xg[nxg];
//                         double dry = gdm->yg[j] - gdm->yg[jrp];
//                         double drz = gdm->zg[0] - gdm->zg[krp];
//                         double dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_z[i][j][0] += dudx1 / dr;

//                         // Upper boundary in z-direction
//                         drz = gdm->zg[nzg + 1] - gdm->zg[krp];
//                         dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_z[i][j][1] += dudx1 / dr;
//                     }
//                 }
//             }
//         }
//     }
    


//     // Contributions from boundaries in y-direction
//     if (sdm->is_y0_boundary) {
//         double dyg1 = sdm->dyg[1];
//         double dyg2 = sdm->dyg[2];
//         double alpha = dyg2 / dyg1;
//         double w1 = gf_coeff / dyg1 * (2.0 + alpha) / (1.0 + alpha);
//         double w2 = gf_coeff / dyg2 / (1.0 + alpha);

//         for (int kp = 0; kp < nzs; ++kp) {
//             for (int ip = 0; ip < nxs; ++ip) {
//                 double dudy0 = -((u[ip][1][kp] - u[ip][0][kp]) * w1 - (u[ip][2][kp] - u[ip][1][kp]) * w2) * sdm->dxm[ip] * sdm->dzm[kp];
//                 int irp = ip + sdm->ista - 1;
//                 int krp = kp + sdm->ksta - 1;
                
//                 // Boundary treatment in x-direction
//                 for (int k = 0; k < nzg; ++k) {
//                     for (int j = 0; j < nyg; ++j) {
//                         // Lower boundary in x-direction x=0
//                         double drx = gdm->xg[0] - gdm->xg[irp];
//                         double dry = gdm->yg[j] - gdm->yg[1];
//                         double drz = gdm->zg[k] - gdm->zg[krp];
//                         double dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_x[j][k][0] += dudy0 / dr;

//                         // Upper boundary in x-direction x=nx+1
//                         drx = gdm->xg[nxg + 1] - gdm->xg[irp];
//                         dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_x[j][k][1] += dudy0 / dr;
//                     }
//                 }

//                 // Boundary treatment in y-direction
//                 for (int k = 0; k < nzg; ++k) {
//                     for (int i = 0; i < nxg; ++i) {
//                         // Lower boundary in y-direction
//                         double drx = gdm->xg[i] - gdm->xg[irp];
//                         double dry = gdm->yg[0] - gdm->yg[1];
//                         double drz = gdm->zg[k] - gdm->zg[krp];
//                         double dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_y[i][k][0] += dudy0 / dr;

//                         // Upper boundary in y-direction
//                         dry = gdm->yg[nyg + 1] - gdm->yg[1];
//                         dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_y[i][k][1] += dudy0 / dr;
//                     }
//                 }

//                 // Boundary treatment in z-direction
//                 for (int j = 0; j < nyg; ++j) {
//                     for (int i = 0; i < nxg; ++i) {
//                         // Lower boundary in z-direction
//                         double drx = gdm->xg[i] - gdm->xg[irp];
//                         double dry = gdm->yg[j] - gdm->yg[1];
//                         double drz = gdm->zg[0] - gdm->zg[krp];
//                         double dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_z[i][j][0] += dudy0 / dr;

//                         // Upper boundary in z-direction
//                         drz = gdm->zg[nzg + 1] - gdm->zg[krp];
//                         dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_z[i][j][1] += dudy0 / dr;
//                     }
//                 }
//             }
//         }
//     }

//     // Contributions from boundaries in y-direction
//     if (sdm->is_y1_boundary) {
//         double dyg1 = sdm->dyg[nxs-2];
//         double dyg2 = sdm->dyg[nxs-1];
//         double alpha = dyg1 / dyg2;
//         double w1 = gf_coeff / dyg2 * (2.0 + alpha) / (1.0 + alpha);
//         double w2 = gf_coeff / dyg1 / (1.0 + alpha);

//         for (int kp = 0; kp < nzs; ++kp) {
//             for (int ip = 0; ip < nxs; ++ip) {
//                 int irp = ip + sdm->ista - 1;
//                 int krp = kp + sdm->ksta - 1;
//                 double dudy1 = -((u[ip][nys-1][kp] - u[ip][nys-2][kp]) * w1 - (u[ip][nys-2][kp] - u[ip][nys-3][kp]) * w2) * sdm->dxm[ip] * sdm->dzm[kp];

//                 // Boundary treatment in x-direction
//                 for (int k = 0; k < nzg; ++k) {
//                     for (int j = 0; j < nyg; ++j) {
//                         // Lower boundary in x-direction x=0
//                         double drx = gdm->xg[0] - gdm->xg[irp];
//                         double dry = gdm->yg[j] - gdm->yg[nyg];
//                         double drz = gdm->zg[k] - gdm->zg[krp];
//                         double dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_x[j][k][0] += dudy1 / dr;

//                         // Upper boundary in x-direction x=nx+1
//                         drx = gdm->xg[nxg + 1] - gdm->xg[irp];
//                         dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_x[j][k][1] += dudy1 / dr;
//                     }
//                 }

//                 // Boundary treatment in y-direction
//                 for (int k = 0; k < nzg; ++k) {
//                     for (int i = 0; i < nxg; ++i) {
//                         // Lower boundary in y-direction
//                         double drx = gdm->xg[i] - gdm->xg[irp];
//                         double dry = gdm->yg[0] - gdm->yg[nyg];
//                         double drz = gdm->zg[k] - gdm->zg[krp];
//                         double dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_y[i][k][0] += dudy1 / dr;

//                         // Upper boundary in y-direction
//                         dry = gdm->yg[nyg + 1] - gdm->yg[nyg];
//                         dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_y[i][k][1] += dudy1 / dr;
//                     }
//                 }

//                 // Boundary treatment in z-direction
//                 for (int j = 0; j < nyg; ++j) {
//                     for (int i = 0; i < nxg; ++i) {
//                         // Lower boundary in z-direction
//                         double drx = gdm->xg[i] - gdm->xg[irp];
//                         double dry = gdm->yg[j] - gdm->yg[nyg];
//                         double drz = gdm->zg[0] - gdm->zg[krp];
//                         double dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_z[i][j][0] += dudy1 / dr;

//                         // Upper boundary in z-direction
//                         drz = gdm->zg[nzg + 1] - gdm->zg[krp];
//                         dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_z[i][j][1] += dudy1 / dr;
//                     }
//                 }
//             }
//         }
//     }


//     // Contributions from boundaries in z-direction
//     if (sdm->is_z0_boundary) {
//         double dzg1 = sdm->dzg[1];
//         double dzg2 = sdm->dzg[2];
//         double alpha = dzg2 / dzg1;
//         double w1 = gf_coeff / dzg1 * (2.0 + alpha) / (1.0 + alpha);
//         double w2 = gf_coeff / dzg2 / (1.0 + alpha);

//         for (int jp = 0; jp < nys; ++jp) {
//             for (int ip = 0; ip < nxs; ++ip) {
//                 int irp = ip + sdm->ista - 1;
//                 int jrp = jp + sdm->jsta - 1;
                
//                 double dudz0 = -((u[ip][jp][1] - u[ip][jp][0]) * w1 - (u[ip][jp][2] - u[ip][jp][1]) * w2) * sdm->dxm[ip] * sdm->dym[jp];

//                 // Boundary treatment in x-direction
//                 for (int k = 0; k < nzg; ++k) {
//                     for (int j = 0; j < nyg; ++j) {
//                         // Lower boundary in x-direction x=0
//                         double drx = gdm->xg[0] - gdm->xg[irp];
//                         double dry = gdm->yg[j] - gdm->yg[jrp];
//                         double drz = gdm->zg[k] - gdm->zg[1];
//                         double dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_x[j][k][0] += dudz0 / dr;

//                         // Upper boundary in x-direction x=nx+1
//                         drx = gdm->xg[nxg + 1] - gdm->xg[irp];
//                         dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_x[j][k][1] += dudz0 / dr;
//                     }
//                 }

//                 // Boundary treatment in y-direction
//                 for (int k = 0; k < nzg; ++k) {
//                     for (int i = 0; i < nxg; ++i) {
//                         // Lower boundary in y-direction
//                         double drx = gdm->xg[i] - gdm->xg[irp];
//                         double dry = gdm->yg[0] - gdm->yg[jrp];
//                         double drz = gdm->zg[k] - gdm->zg[1];
//                         double dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_y[i][k][0] += dudz0 / dr;

//                         // Upper boundary in y-direction
//                         dry = gdm->yg[nyg + 1] - gdm->yg[jrp];
//                         dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_y[i][k][1] += dudz0 / dr;
//                     }
//                 }

//                 // Boundary treatment in z-direction
//                 for (int j = 0; j < nyg; ++j) {
//                     for (int i = 0; i < nxg; ++i) {
//                         // Lower boundary in z-direction
//                         double drx = gdm->xg[i] - gdm->xg[irp];
//                         double dry = gdm->yg[j] - gdm->yg[jrp];
//                         double drz = gdm->zg[0] - gdm->zg[1];
//                         double dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_z[i][j][0] += dudz0 / dr;

//                         // Upper boundary in z-direction
//                         drz = gdm->zg[nzg + 1] - gdm->zg[1];
//                         dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_z[i][j][1] += dudz0 / dr;
//                     }
//                 }
//             }
//         }
//     }

//     // Contributions from boundaries in z-direction
//     if (sdm->is_z1_boundary) {
//         double dzg1 = sdm->dzg[nzs-2];
//         double dzg2 = sdm->dzg[nzs-1];
//         double alpha = dzg1 / dzg2;
//         double w1 = gf_coeff / dzg2 * (2.0 + alpha) / (1.0 + alpha);
//         double w2 = gf_coeff / dzg1 / (1.0 + alpha);

//         for (int jp = 0; jp < nys; ++jp) {
//             for (int ip = 0; ip < nxs; ++ip) {
//                 int irp = ip + sdm->ista - 1;
//                 int jrp = jp + sdm->jsta - 1;
                
//                 double dudz1 = -((u[ip][jp][nzs-1] - u[ip][jp][nzs-2]) * w1 - (u[ip][jp][nzs-2] - u[ip][jp][nzs-3]) * w2) * sdm->dxm[ip] * sdm->dym[jp];

//                 // Boundary treatment in x-direction
//                 for (int k = 0; k < nzg; ++k) {
//                     for (int j = 0; j < nyg; ++j) {
//                         // Lower boundary in x-direction x=0
//                         double drx = gdm->xg[0] - gdm->xg[irp];
//                         double dry = gdm->yg[j] - gdm->yg[jrp];
//                         double drz = gdm->zg[k] - gdm->zg[nzg];
//                         double dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_x[j][k][0] += dudz1 / dr;

//                         // Upper boundary in x-direction x=nx+1
//                         drx = gdm->xg[nxg + 1] - gdm->xg[irp];
//                         dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_x[j][k][1] += dudz1 / dr;
//                     }
//                 }

//                 // Boundary treatment in y-direction
//                 for (int k = 0; k < nzg; ++k) {
//                     for (int i = 0; i < nxg; ++i) {
//                         // Lower boundary in y-direction
//                         double drx = gdm->xg[i] - gdm->xg[irp];
//                         double dry = gdm->yg[0] - gdm->yg[jrp];
//                         double drz = gdm->zg[k] - gdm->zg[nzg];
//                         double dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_y[i][k][0] += dudz1 / dr;

//                         // Upper boundary in y-direction
//                         dry = gdm->yg[nyg + 1] - gdm->yg[jrp];
//                         dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_y[i][k][1] += dudz1 / dr;
//                     }
//                 }

//                 // Boundary treatment in z-direction
//                 for (int j = 0; j < nyg; ++j) {
//                     for (int i = 0; i < nxg; ++i) {
//                         // Lower boundary in z-direction
//                         double drx = gdm->xg[i] - gdm->xg[irp];
//                         double dry = gdm->yg[j] - gdm->yg[jrp];
//                         double drz = gdm->zg[0] - gdm->zg[nzg];
//                         double dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_z[i][j][0] += dudz1 / dr;

//                         // Upper boundary in z-direction
//                         drz = gdm->zg[nzg + 1] - gdm->zg[nzg];
//                         dr = sqrt(drx*drx + dry*dry + drz*drz);
//                         bval_z[i][j][1] += dudz1 / dr;
//                     }
//                 }
//             }
//         }
//     }


// #ifdef MPI_INPLACE
//     if (comm_boundary->mpi_comm != MPI_COMM_NULL) {
//         MPI_Allreduce(MPI_IN_PLACE, bval_x, bval_x_size, MPI_DOUBLE, MPI_SUM, comm_boundary->mpi_comm);
//         MPI_Allreduce(MPI_IN_PLACE, bval_y, bval_y_size, MPI_DOUBLE, MPI_SUM, comm_boundary->mpi_comm);
//         MPI_Allreduce(MPI_IN_PLACE, bval_z, bval_z_size, MPI_DOUBLE, MPI_SUM, comm_boundary->mpi_comm);
//     }

//     if (sdm->is_x0_boundary) {
//         for (int k = 0; k < nzs; ++k)
//             for (int j = 0; j < nys; ++j)
//                 u[0][j][k] = -bval_x[j+sdm->jsta-1][k+sdm->ksta-1][0];
//     }
//     if (sdm->is_x1_boundary) {
//         for (int k = 0; k < nzs; ++k)
//             for (int j = 0; j < nys; ++j)
//                 u[nxs+1][j][k] = -bval_x[j+sdm->jsta-1][k+sdm->ksta-1][1];
//     }
//     if (sdm->is_y0_boundary) {
//         for (int k = 0; k < nzs; ++k)
//             for (int i = 0; i < nxs; ++i)
//                 u[i][0][k] = -bval_y[i+sdm->ista-1][k+sdm->ksta-1][0];
//     }
//     if (sdm->is_y1_boundary) {
//         for (int k = 0; k < nzs; ++k)
//             for (int i = 0; i < nxs; ++i)
//                 u[i][nys+1][k] = -bval_y[i+sdm->ista-1][k+sdm->ksta-1][1];
//     }
//     if (sdm->is_z0_boundary) {
//         for (int j = 0; j < nys; ++j)
//             for (int i = 0; i < nxs; ++i)
//                 u[i][j][0] = -bval_z[i+sdm->ista-1][j+sdm->jsta-1][0];
//     }
//     if (sdm->is_z1_boundary) {
//         for (int j = 0; j < nys; ++j)
//             for (int i = 0; i < nxs; ++i)
//                 u[i][j][nzs+1] = -bval_z[i+sdm->ista-1][j+sdm->jsta-1][1];
//     }
//     if (myrank == 0) printf("[Boundary calculation] With MPI_IN_PLACE\n");
// #else
//     if (comm_boundary->mpi_comm != MPI_COMM_NULL) {
//         MPI_Allreduce(bval_x, bval_xg, (nyg)*(nzg)*2, MPI_DOUBLE, MPI_SUM, comm_boundary->mpi_comm);
//         MPI_Allreduce(bval_y, bval_yg, (nxg)*(nzg)*2, MPI_DOUBLE, MPI_SUM, comm_boundary->mpi_comm);
//         MPI_Allreduce(bval_z, bval_zg, (nxg)*(nyg)*2, MPI_DOUBLE, MPI_SUM, comm_boundary->mpi_comm);
//     }

//     if (sdm->is_x0_boundary) {
//         for (int k = 0; k < nzs; ++k)
//             for (int j = 0; j < nys; ++j)
//                 u[0][j][k] = -bval_xg[j+sdm->jsta-1][k+sdm->ksta-1][0];
//     }
//     if (sdm->is_x1_boundary) {
//         for (int k = 0; k < nzs; ++k)
//             for (int j = 0; j < nys; ++j)
//                 u[nxs+1][j][k] = -bval_xg[j+sdm->jsta-1][k+sdm->ksta-1][1];
//     }
//     if (sdm->is_y0_boundary) {
//         for (int k = 0; k < nzs; ++k)
//             for (int i = 0; i < nxs; ++i)
//                 u[i][0][k] = -bval_yg[i+sdm->ista-1][k+sdm->ksta-1][0];
//     }
//     if (sdm->is_y1_boundary) {
//         for (int k = 0; k < nzs; ++k)
//             for (int i = 0; i < nxs; ++i)
//                 u[i][nys+1][k] = -bval_yg[i+sdm->ista-1][k+sdm->ksta-1][1];
//     }
//     if (sdm->is_z0_boundary) {
//         for (int j = 0; j < nys; ++j)
//             for (int i = 0; i < nxs; ++i)
//                 u[i][j][0] = -bval_zg[i+sdm->ista-1][j+sdm->jsta-1][0];
//     }
//     if (sdm->is_z1_boundary) {
//         for (int j = 0; j < nys; ++j)
//             for (int i = 0; i < nxs; ++i)
//                 u[i][j][nzs+1] = -bval_zg[i+sdm->ista-1][j+sdm->jsta-1][1];
//     }
//     if (myrank == 0) printf("[Boundary calculation] Without MPI_IN_PLACE\n");

//     free(bval_xg);
//     free(bval_yg);
//     free(bval_zg);
// #endif

// free(bval_x);
// free(bval_y);
// free(bval_z);

// }

// void geometry_boundary_normal_derivative(
//     double ***dudx, double ***dudy, double ***dudz,
//     double ***u, subdomain *sdm)
// {
//     int i, j, k;
//     int nx = sdm->nx;
//     int ny = sdm->ny;
//     int nz = sdm->nz;

//     double dxg1, dxg2, dyg1, dyg2, dzg1, dzg2, alpha;

//     if (sdm->is_x0_boundary) {
//         dxg1 = sdm->dxg[2];
//         dxg2 = sdm->dxg[3];
//         alpha = dxg2 / dxg1;

//         for (k = 1; k <= nz; ++k)
//             for (j = 1; j <= ny; ++j)
//                 dudx[0][j][k] = -(
//                     (u[2][j][k] - u[1][j][k]) / dxg1 * (2.0 + alpha) / (1.0 + alpha)
//                     - (u[3][j][k] - u[2][j][k]) / dxg2 / (1.0 + alpha)
//                 );
//     }

//     if (sdm->is_x1_boundary) {
//         dxg1 = sdm->dxg[nx - 1];
//         dxg2 = sdm->dxg[nx];
//         alpha = dxg1 / dxg2;

//         for (k = 1; k <= nz; ++k)
//             for (j = 1; j <= ny; ++j)
//                 dudx[1][j][k] = (
//                     (u[nx][j][k] - u[nx - 1][j][k]) / dxg2 * (2.0 + alpha) / (1.0 + alpha)
//                     - (u[nx - 1][j][k] - u[nx - 2][j][k]) / dxg1 / (1.0 + alpha)
//                 );
//     }

//     if (sdm->is_y0_boundary) {
//         dyg1 = sdm->dyg[2];
//         dyg2 = sdm->dyg[3];
//         alpha = dyg2 / dyg1;

//         for (k = 1; k <= nz; ++k)
//             for (i = 1; i <= nx; ++i)
//                 dudy[i][0][k] = -(
//                     (u[i][2][k] - u[i][1][k]) / dyg1 * (2.0 + alpha) / (1.0 + alpha)
//                     - (u[i][3][k] - u[i][2][k]) / dyg2 / (1.0 + alpha)
//                 );
//     }

//     if (sdm->is_y1_boundary) {
//         dyg1 = sdm->dyg[ny - 1];
//         dyg2 = sdm->dyg[ny];
//         alpha = dyg1 / dyg2;

//         for (k = 1; k <= nz; ++k)
//             for (i = 1; i <= nx; ++i)
//                 dudy[i][1][k] = (
//                     (u[i][ny][k] - u[i][ny - 1][k]) / dyg2 * (2.0 + alpha) / (1.0 + alpha)
//                     - (u[i][ny - 1][k] - u[i][ny - 2][k]) / dyg1 / (1.0 + alpha)
//                 );
//     }

//     if (sdm->is_z0_boundary) {
//         dzg1 = sdm->dzg[2];
//         dzg2 = sdm->dzg[3];
//         alpha = dzg2 / dzg1;

//         for (j = 1; j <= ny; ++j)
//             for (i = 1; i <= nx; ++i)
//                 dudz[i][j][0] = -(
//                     (u[i][j][2] - u[i][j][1]) / dzg1 * (2.0 + alpha) / (1.0 + alpha)
//                     - (u[i][j][3] - u[i][j][2]) / dzg2 / (1.0 + alpha)
//                 );
//     }

//     if (sdm->is_z1_boundary) {
//         dzg1 = sdm->dzg[nz - 1];
//         dzg2 = sdm->dzg[nz];
//         alpha = dzg1 / dzg2;

//         for (j = 1; j <= ny; ++j)
//             for (i = 1; i <= nx; ++i)
//                 dudz[i][j][1] = (
//                     (u[i][j][nz] - u[i][j][nz - 1]) / dzg2 * (2.0 + alpha) / (1.0 + alpha)
//                     - (u[i][j][nz - 1] - u[i][j][nz - 2]) / dzg1 / (1.0 + alpha)
//                 );
//     }
// }

// void geometry_boundary_surface_integral(
//     double ***u, double ***dudx, double ***dudy, double ***dudz,
//     subdomain *sdm, domain *gdm) {

//     int myrank;
//     MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

//     int nxg = gdm->nx;
//     int nyg = gdm->ny;
//     int nzg = gdm->nz;

//     int nxs = sdm->nx;
//     int nys = sdm->ny;
//     int nzs = sdm->nz;

//     // Allocate and initialize bval arrays (3D arrays with 3rd dim of size 2)
//     double ***bval_x = calloc(nyg, sizeof(double**));
//     double ***bval_y = calloc(nxg, sizeof(double**));
//     double ***bval_z = calloc(nxg, sizeof(double**));

//     double ***bval_xg = calloc(nyg, sizeof(double**));
//     double ***bval_yg = calloc(nxg, sizeof(double**));
//     double ***bval_zg = calloc(nxg, sizeof(double**));

//     for (int j = 0; j < nyg; ++j) {
//         bval_x[j] = calloc(nzg, sizeof(double*));
//         bval_xg[j] = calloc(nzg, sizeof(double*));
//         for (int k = 0; k < nzg; ++k) {
//             bval_x[j][k] = calloc(2, sizeof(double));
//             bval_xg[j][k] = calloc(2, sizeof(double));
//         }
//     }

//     for (int i = 0; i < nxg; ++i) {
//         bval_y[i] = calloc(nzg, sizeof(double*));
//         bval_yg[i] = calloc(nzg, sizeof(double*));
//         bval_z[i] = calloc(nyg, sizeof(double*));
//         bval_zg[i] = calloc(nyg, sizeof(double*));
//         for (int k = 0; k < nzg; ++k) {
//             bval_y[i][k] = calloc(2, sizeof(double));
//             bval_yg[i][k] = calloc(2, sizeof(double));
//         }
//         for (int j = 0; j < nyg; ++j) {
//             bval_z[i][j] = calloc(2, sizeof(double));
//             bval_zg[i][j] = calloc(2, sizeof(double));
//         }
//     }

 
//     int ierr;
//     // Flatten arrays for MPI_Allreduce
//     int size_x = nyg * nzg * 2;
//     int size_y = nxg * nzg * 2;
//     int size_z = nxg * nyg * 2;
//     double *flat_x = malloc(sizeof(double) * size_x);
//     double *flat_xg = calloc(size_x, sizeof(double));
//     double *flat_y = malloc(sizeof(double) * size_y);
//     double *flat_yg = calloc(size_y, sizeof(double));
//     double *flat_z = malloc(sizeof(double) * size_z);
//     double *flat_zg = calloc(size_z, sizeof(double));

//     // Flatten
//     int idx = 0;
//     for (int j = 0; j < nyg; ++j)
//         for (int k = 0; k < nzg; ++k)
//             for (int b = 0; b < 2; ++b)
//                 flat_x[idx++] = bval_x[j][k][b];

//     idx = 0;
//     for (int i = 0; i < nxg; ++i)
//         for (int k = 0; k < nzg; ++k)
//             for (int b = 0; b < 2; ++b)
//                 flat_y[idx++] = bval_y[i][k][b];

//     idx = 0;
//     for (int i = 0; i < nxg; ++i)
//         for (int j = 0; j < nyg; ++j)
//             for (int b = 0; b < 2; ++b)
//                 flat_z[idx++] = bval_z[i][j][b];

//     // Reduce
//     MPI_Allreduce(flat_x, flat_xg, size_x, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//     MPI_Allreduce(flat_y, flat_yg, size_y, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//     MPI_Allreduce(flat_z, flat_zg, size_z, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

//     // Scatter back into 3D arrays
//     idx = 0;
//     for (int j = 0; j < nyg; ++j)
//         for (int k = 0; k < nzg; ++k)
//             for (int b = 0; b < 2; ++b)
//                 bval_xg[j][k][b] = flat_xg[idx++];

//     idx = 0;
//     for (int i = 0; i < nxg; ++i)
//         for (int k = 0; k < nzg; ++k)
//             for (int b = 0; b < 2; ++b)
//                 bval_yg[i][k][b] = flat_yg[idx++];

//     idx = 0;
//     for (int i = 0; i < nxg; ++i)
//         for (int j = 0; j < nyg; ++j)
//             for (int b = 0; b < 2; ++b)
//                 bval_zg[i][j][b] = flat_zg[idx++];

//     // Write to file (one file per rank)
//     char myfilename[256];
//     snprintf(myfilename, sizeof(myfilename), "./result/bval.%03d", myrank);
//     FILE *fout = fopen(myfilename, "w");

//     fprintf(fout, "bval_xg\n");
//     for (int k = 0; k < nzg; ++k)
//         for (int j = 0; j < nyg; ++j)
//             fprintf(fout, "%9.5f %9.5f %15.7e %15.7e\n", gdm->yg[j], gdm->zg[k], bval_xg[j][k][0], bval_xg[j][k][1]);

//     fprintf(fout, "bval_yg\n");
//     for (int k = 0; k < nzg; ++k)
//         for (int i = 0; i < nxg; ++i)
//             fprintf(fout, "%9.5f %9.5f %15.7e %15.7e\n", gdm->xg[i], gdm->zg[k], bval_yg[i][k][0], bval_yg[i][k][1]);

//     fprintf(fout, "bval_zg\n");
//     for (int j = 0; j < nyg; ++j)
//         for (int i = 0; i < nxg; ++i)
//             fprintf(fout, "%9.5f %9.5f %15.7e %15.7e\n", gdm->xg[i], gdm->yg[j], bval_zg[i][j][0], bval_zg[i][j][1]);

//     fclose(fout);

//     if (sdm->is_x0_boundary)
//     {
//         for (int k = 0; k < nzs; k++)
//         {
//             for (int j = 0; j < nys; j++)
//             {
//                 u[0][j][k] =  -bval_xg[j+sdm->jsta-1][k+sdm->ksta-1][0];
//             }
//         }
//     }
//     if (sdm->is_x1_boundary)
//     {
//         for (int k = 0; k < nzs; k++)
//         {
//             for (int j = 0; j < nys; j++)
//             {
//                 u[nxs+1][j][k] =  -bval_xg[j+sdm->jsta-1][k+sdm->ksta-1][1];
//             }
//         }
//     }

//     if (sdm->is_y0_boundary)
//     {
//         for (int k = 0; k < nzs; k++)
//         {
//             for (int i = 0; i < nxs; i++)
//             {
//                 u[i][0][k] =  -bval_yg[i+sdm->ista-1][k+sdm->ksta-1][0];
//             }
//         }
//     }
//     if (sdm->is_y1_boundary)
//     {
//         for (int k = 0; k < nzs; k++)
//         {
//             for (int i = 0; i < nxs; i++)
//             {
//                 u[i][nys+1][k] =  -bval_yg[i+sdm->ista-1][k+sdm->ksta-1][1];
//             }
//         }
//     }

//     if (sdm->is_z0_boundary)
//     {
//         for (int j = 0; j < nys; j++)
//         {
//             for (int i = 0; i < nxs; i++)
//             {
//                 u[i][j][0] =  -bval_zg[i+sdm->ista-1][j+sdm->jsta-1][0];
//             }
//         }
//     }
//     if (sdm->is_z1_boundary)
//     {
//         for (int j = 0; j < nys; j++)
//         {
//             for (int i = 0; i < nxs; i++)
//             {
//                 u[i][j][nzs+1] =  -bval_zg[i+sdm->ista-1][j+sdm->jsta-1][1];
//             }
//         }
//     }
    

//     // Free memory
//     free(flat_x); free(flat_y); free(flat_z);
//     free(flat_xg); free(flat_yg); free(flat_zg);

    
// }
