#include "geometry.h"
#include <cuda_runtime.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#ifndef IDX
#define IDX(i,j,k,ni,nj) ((i)*(ni) + (j)*(nj) + (k))
#endif

#define CUDA_CHECK(call) do {                                      \
    cudaError_t err__ = (call);                                    \
    if (err__ != cudaSuccess) {                                    \
        fprintf(stderr, "CUDA error at %s:%d: %s\n",               \
                __FILE__, __LINE__, cudaGetErrorString(err__));    \
        exit(EXIT_FAILURE);                                        \
    }                                                             \
} while (0)

#include <limits.h>


static double *halo_send_low = NULL;
static double *halo_send_high = NULL;
static double *halo_recv_low = NULL;
static double *halo_recv_high = NULL;
static size_t halo_capacity = 0;

static void geometry_halo_buffers_release_gpu(void)
{
    if (halo_send_low)  CUDA_CHECK(cudaFree(halo_send_low));
    if (halo_send_high) CUDA_CHECK(cudaFree(halo_send_high));
    if (halo_recv_low)  CUDA_CHECK(cudaFree(halo_recv_low));
    if (halo_recv_high) CUDA_CHECK(cudaFree(halo_recv_high));
    halo_send_low = halo_send_high = NULL;
    halo_recv_low = halo_recv_high = NULL;
    halo_capacity = 0;
}

static void geometry_halo_buffers_ensure_gpu(size_t count)
{
    if (count <= halo_capacity) return;

    geometry_halo_buffers_release_gpu();
    CUDA_CHECK(cudaMalloc((void **)&halo_send_low,  count * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void **)&halo_send_high, count * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void **)&halo_recv_low,  count * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void **)&halo_recv_high, count * sizeof(double)));
    halo_capacity = count;
}

__global__ static void geometry_pack_face_kernel(
    const double *u, double *low, double *high,
    int nx, int ny, int nz, int direction, size_t count)
{
    size_t tid = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= count) return;

    int i, j, k;
    int ni = (ny + 2) * (nz + 2);
    int nj = nz + 2;
    size_t low_idx, high_idx;

    if (direction == 0) {
        j = (int)(tid / nz) + 1;
        k = (int)(tid % nz) + 1;
        low_idx  = (size_t)1  * ni + (size_t)j * nj + k;
        high_idx = (size_t)nx * ni + (size_t)j * nj + k;
    } else if (direction == 1) {
        i = (int)(tid / nz) + 1;
        k = (int)(tid % nz) + 1;
        low_idx  = (size_t)i * ni + (size_t)1  * nj + k;
        high_idx = (size_t)i * ni + (size_t)ny * nj + k;
    } else {
        i = (int)(tid / ny) + 1;
        j = (int)(tid % ny) + 1;
        low_idx  = (size_t)i * ni + (size_t)j * nj + 1;
        high_idx = (size_t)i * ni + (size_t)j * nj + nz;
    }

    low[tid] = u[low_idx];
    high[tid] = u[high_idx];
}

__global__ static void geometry_unpack_face_kernel(
    double *u, const double *low, const double *high,
    int nx, int ny, int nz, int direction, size_t count,
    int unpack_low, int unpack_high)
{
    size_t tid = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= count) return;

    int i, j, k;
    int ni = (ny + 2) * (nz + 2);
    int nj = nz + 2;
    size_t low_idx, high_idx;

    if (direction == 0) {
        j = (int)(tid / nz) + 1;
        k = (int)(tid % nz) + 1;
        low_idx  = (size_t)0        * ni + (size_t)j * nj + k;
        high_idx = (size_t)(nx + 1) * ni + (size_t)j * nj + k;
    } else if (direction == 1) {
        i = (int)(tid / nz) + 1;
        k = (int)(tid % nz) + 1;
        low_idx  = (size_t)i * ni + (size_t)0        * nj + k;
        high_idx = (size_t)i * ni + (size_t)(ny + 1) * nj + k;
    } else {
        i = (int)(tid / ny) + 1;
        j = (int)(tid % ny) + 1;
        low_idx  = (size_t)i * ni + (size_t)j * nj;
        high_idx = (size_t)i * ni + (size_t)j * nj + (nz + 1);
    }

    if (unpack_low)  u[low_idx] = low[tid];
    if (unpack_high) u[high_idx] = high[tid];
}

static void geometry_exchange_one_direction_gpu(
    double *u, subdomain *sdm, int direction,
    const cart_comm_1d *comm, int tag_low_to_high, int tag_high_to_low)
{
    size_t count = direction == 0 ? (size_t)sdm->ny * sdm->nz
                 : direction == 1 ? (size_t)sdm->nx * sdm->nz
                                  : (size_t)sdm->nx * sdm->ny;
    if (count == 0) return;
    if (count > (size_t)INT_MAX) {
        fprintf(stderr, "halo MPI count exceeds INT_MAX\n");
        MPI_Abort(comm->mpi_comm, EXIT_FAILURE);
    }

    geometry_halo_buffers_ensure_gpu(count);

    const int threads = 256;
    const int blocks = (int)((count + threads - 1) / threads);
    geometry_pack_face_kernel<<<blocks, threads>>>(
        u, halo_send_low, halo_send_high,
        sdm->nx, sdm->ny, sdm->nz, direction, count);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    MPI_Request request[4];
    int nrequest = 0;

    if (comm->east_rank != MPI_PROC_NULL) {
        MPI_Isend(halo_send_high, (int)count, MPI_DOUBLE,
                  comm->east_rank, tag_low_to_high, comm->mpi_comm,
                  &request[nrequest++]);
        MPI_Irecv(halo_recv_high, (int)count, MPI_DOUBLE,
                  comm->east_rank, tag_high_to_low, comm->mpi_comm,
                  &request[nrequest++]);
    }
    if (comm->west_rank != MPI_PROC_NULL) {
        MPI_Isend(halo_send_low, (int)count, MPI_DOUBLE,
                  comm->west_rank, tag_high_to_low, comm->mpi_comm,
                  &request[nrequest++]);
        MPI_Irecv(halo_recv_low, (int)count, MPI_DOUBLE,
                  comm->west_rank, tag_low_to_high, comm->mpi_comm,
                  &request[nrequest++]);
    }
    if (nrequest) MPI_Waitall(nrequest, request, MPI_STATUSES_IGNORE);

    const int unpack_low = comm->west_rank != MPI_PROC_NULL;
    const int unpack_high = comm->east_rank != MPI_PROC_NULL;
    if (unpack_low || unpack_high) {
        geometry_unpack_face_kernel<<<blocks, threads>>>(
            u, halo_recv_low, halo_recv_high,
            sdm->nx, sdm->ny, sdm->nz, direction, count,
            unpack_low, unpack_high);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
    }
}



// __global__
static void create_1d_geometry_host(double *dm, double *dg, double *xg, int n, double origin, double length, double stretching)
{
    int i;
    if (stretching == 1.0) {
        dm[0] = 0.0;
        dg[0] = 0.0;
        xg[0] = origin;

        for (i = 1; i <= n; i++) {
            dm[i] = length / (double)n;
            dg[i] = 0.5 * (dm[i] + dm[i-1]);
            xg[i] = xg[i-1] + dg[i];
        }

        dm[n+1] = 0.0;
        dg[n+1] = 0.5 * (dm[n+1] + dm[n]);
        xg[n+1] = xg[n] + dg[n+1];
    } else {
        dm[0] = 0.0;

        dm[1] = 0.5 * length * (stretching - 1.0) / (pow(stretching, n/2) - 1.0);

        for (i = 2; i <= n/2; i++) {
            dm[i] = dm[i-1] * stretching;
            dm[n - i + 1] = dm[i];
        }

        dm[n] = dm[1];
        dm[n+1] = 0.0;

        dg[0] = 0.0;
        xg[0] = origin;

        for (i = 1; i <= n; i++) {
            dg[i] = 0.5 * (dm[i] + dm[i-1]);
            xg[i] = xg[i-1] + dg[i];
        }

        dg[n+1] = 0.5 * (dm[n+1] + dm[n]);
        xg[n+1] = xg[n] + dg[n+1];
    }
}

void geometry_domain_create_gpu(domain *gdm,
                            int nx, int ny, int nz,
                            double ox, double oy, double oz,
                            double lx, double ly, double lz,
                            double ax, double ay, double az,
                            int period[3])
{
    gdm->nx = nx;
    gdm->ny = ny;
    gdm->nz = nz;

    gdm->ox = ox;
    gdm->oy = oy;
    gdm->oz = oz;

    gdm->lx = lx;
    gdm->ly = ly;
    gdm->lz = lz;

    for (int d = 0; d < 3; d++) {
        gdm->is_periodic[d] = period[d];
    }

    double *h_dxm = (double*)calloc(nx + 2, sizeof(double));
    double *h_dym = (double*)calloc(ny + 2, sizeof(double));
    double *h_dzm = (double*)calloc(nz + 2, sizeof(double));

    double *h_dxg = (double*)calloc(nx + 2, sizeof(double));
    double *h_dyg = (double*)calloc(ny + 2, sizeof(double));
    double *h_dzg = (double*)calloc(nz + 2, sizeof(double));

    double *h_xg  = (double*)calloc(nx + 2, sizeof(double));
    double *h_yg  = (double*)calloc(ny + 2, sizeof(double));
    double *h_zg  = (double*)calloc(nz + 2, sizeof(double));

    if (!h_dxm || !h_dym || !h_dzm ||
        !h_dxg || !h_dyg || !h_dzg ||
        !h_xg  || !h_yg  || !h_zg) {
        fprintf(stderr, "Error: host memory allocation failed in geometry_domain_create\n");
        exit(EXIT_FAILURE);
    }

    create_1d_geometry_host(h_dxm, h_dxg, h_xg, nx, ox, lx, ax);
    create_1d_geometry_host(h_dym, h_dyg, h_yg, ny, oy, ly, ay);
    create_1d_geometry_host(h_dzm, h_dzg, h_zg, nz, oz, lz, az);

    CUDA_CHECK(cudaMalloc((void**)&gdm->dxm, (nx + 2) * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&gdm->dym, (ny + 2) * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&gdm->dzm, (nz + 2) * sizeof(double)));

    CUDA_CHECK(cudaMalloc((void**)&gdm->dxg, (nx + 2) * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&gdm->dyg, (ny + 2) * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&gdm->dzg, (nz + 2) * sizeof(double)));

    CUDA_CHECK(cudaMalloc((void**)&gdm->xg,  (nx + 2) * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&gdm->yg,  (ny + 2) * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&gdm->zg,  (nz + 2) * sizeof(double)));

    CUDA_CHECK(cudaMemcpy(gdm->dxm, h_dxm, (nx + 2) * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(gdm->dym, h_dym, (ny + 2) * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(gdm->dzm, h_dzm, (nz + 2) * sizeof(double), cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMemcpy(gdm->dxg, h_dxg, (nx + 2) * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(gdm->dyg, h_dyg, (ny + 2) * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(gdm->dzg, h_dzg, (nz + 2) * sizeof(double), cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMemcpy(gdm->xg, h_xg, (nx + 2) * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(gdm->yg, h_yg, (ny + 2) * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(gdm->zg, h_zg, (nz + 2) * sizeof(double), cudaMemcpyHostToDevice));

    free(h_dxm); free(h_dym); free(h_dzm);
    free(h_dxg); free(h_dyg); free(h_dzg);
    free(h_xg);  free(h_yg);  free(h_zg);
}

void geometry_domain_destroy_gpu(domain *gdm)
{
    if (gdm->dxm) CUDA_CHECK(cudaFree(gdm->dxm));
    if (gdm->dym) CUDA_CHECK(cudaFree(gdm->dym));
    if (gdm->dzm) CUDA_CHECK(cudaFree(gdm->dzm));

    if (gdm->dxg) CUDA_CHECK(cudaFree(gdm->dxg));
    if (gdm->dyg) CUDA_CHECK(cudaFree(gdm->dyg));
    if (gdm->dzg) CUDA_CHECK(cudaFree(gdm->dzg));

    if (gdm->xg)  CUDA_CHECK(cudaFree(gdm->xg));
    if (gdm->yg)  CUDA_CHECK(cudaFree(gdm->yg));
    if (gdm->zg)  CUDA_CHECK(cudaFree(gdm->zg));

    gdm->dxm = gdm->dym = gdm->dzm = NULL;
    gdm->dxg = gdm->dyg = gdm->dzg = NULL;
    gdm->xg  = gdm->yg  = gdm->zg  = NULL;
    geometry_halo_buffers_release_gpu();
}

void geometry_subdomain_create_gpu(subdomain *sdm, const domain *gdm)
{
    for (int d = 0; d < 3; d++) {
        sdm->is_periodic[d]   = gdm->is_periodic[d];
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

    int npx = sdm->nx + 2;
    int npy = sdm->ny + 2;
    int npz = sdm->nz + 2;

    CUDA_CHECK(cudaMalloc((void**)&sdm->dxm, npx * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&sdm->dym, npy * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&sdm->dzm, npz * sizeof(double)));

    CUDA_CHECK(cudaMalloc((void**)&sdm->dxg, npx * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&sdm->dyg, npy * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&sdm->dzg, npz * sizeof(double)));

    CUDA_CHECK(cudaMalloc((void**)&sdm->xg,  npx * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&sdm->yg,  npy * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&sdm->zg,  npz * sizeof(double)));

    CUDA_CHECK(cudaMemcpy(sdm->dxm, gdm->dxm + sdm->ista - 1, npx * sizeof(double), cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaMemcpy(sdm->dxg, gdm->dxg + sdm->ista - 1, npx * sizeof(double), cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaMemcpy(sdm->xg, gdm->xg + sdm->ista - 1, npx * sizeof(double), cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaMemcpy(sdm->dym, gdm->dym + sdm->jsta - 1, npy * sizeof(double), cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaMemcpy(sdm->dyg, gdm->dyg + sdm->jsta - 1, npy * sizeof(double), cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaMemcpy(sdm->yg, gdm->yg + sdm->jsta - 1, npy * sizeof(double), cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaMemcpy(sdm->dzm, gdm->dzm + sdm->ksta - 1, npz * sizeof(double), cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaMemcpy(sdm->dzg, gdm->dzg + sdm->ksta - 1, npz * sizeof(double), cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaMemcpy(sdm->zg, gdm->zg + sdm->ksta - 1, npz * sizeof(double), cudaMemcpyDeviceToDevice));

    double h_xg1, h_dxm1;
    double h_yg1, h_dym1;
    double h_zg1, h_dzm1;

    CUDA_CHECK(cudaMemcpy(&h_xg1,  sdm->xg  + 1, sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(&h_dxm1, sdm->dxm + 1, sizeof(double), cudaMemcpyDeviceToHost));

    CUDA_CHECK(cudaMemcpy(&h_yg1,  sdm->yg  + 1, sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(&h_dym1, sdm->dym + 1, sizeof(double), cudaMemcpyDeviceToHost));

    CUDA_CHECK(cudaMemcpy(&h_zg1,  sdm->zg  + 1, sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(&h_dzm1, sdm->dzm + 1, sizeof(double), cudaMemcpyDeviceToHost));
    
    sdm->ox = h_xg1 - 0.5 * h_dxm1;
    sdm->oy = h_yg1 - 0.5 * h_dym1;
    sdm->oz = h_zg1 - 0.5 * h_dzm1;



    size_t size3d = (size_t)npx * (size_t)npy * (size_t)npz;

    CUDA_CHECK(cudaMalloc((void**)&sdm->x, size3d * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&sdm->b, size3d * sizeof(double)));
    CUDA_CHECK(cudaMalloc((void**)&sdm->r, size3d * sizeof(double)));

    CUDA_CHECK(cudaMemset(sdm->x, 0, size3d * sizeof(double)));
    CUDA_CHECK(cudaMemset(sdm->b, 0, size3d * sizeof(double)));
    CUDA_CHECK(cudaMemset(sdm->r, 0, size3d * sizeof(double)));

    sdm->is_x0_boundary = (comm_1d_x.west_rank == MPI_PROC_NULL);
    sdm->is_x1_boundary = (comm_1d_x.east_rank == MPI_PROC_NULL);

    sdm->is_y0_boundary = (comm_1d_y.west_rank == MPI_PROC_NULL);
    sdm->is_y1_boundary = (comm_1d_y.east_rank == MPI_PROC_NULL);

    sdm->is_z0_boundary = (comm_1d_z.west_rank == MPI_PROC_NULL);
    sdm->is_z1_boundary = (comm_1d_z.east_rank == MPI_PROC_NULL);


}

void geometry_subdomain_destroy_gpu(subdomain *sdm)
{
    if (sdm->dxm) CUDA_CHECK(cudaFree(sdm->dxm));
    if (sdm->dym) CUDA_CHECK(cudaFree(sdm->dym));
    if (sdm->dzm) CUDA_CHECK(cudaFree(sdm->dzm));

    if (sdm->dxg) CUDA_CHECK(cudaFree(sdm->dxg));
    if (sdm->dyg) CUDA_CHECK(cudaFree(sdm->dyg));
    if (sdm->dzg) CUDA_CHECK(cudaFree(sdm->dzg));

    if (sdm->xg) CUDA_CHECK(cudaFree(sdm->xg));
    if (sdm->yg) CUDA_CHECK(cudaFree(sdm->yg));
    if (sdm->zg) CUDA_CHECK(cudaFree(sdm->zg));

    if (sdm->x) CUDA_CHECK(cudaFree(sdm->x));
    if (sdm->b) CUDA_CHECK(cudaFree(sdm->b));
    if (sdm->r) CUDA_CHECK(cudaFree(sdm->r));

    sdm->dxm = sdm->dym = sdm->dzm = NULL;
    sdm->dxg = sdm->dyg = sdm->dzg = NULL;
    sdm->xg  = sdm->yg  = sdm->zg  = NULL;
    sdm->x = sdm->b = sdm->r = NULL;
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

void geometry_halocell_update_selectively_gpu(
    double *u, subdomain *sdm, int *is_serial)
{
    /*
     * Each active direction is packed into contiguous device memory.
     * For the current two-GPU runs normally only one branch is active.
     */
    if (!is_serial[0])
        geometry_exchange_one_direction_gpu(
            u, sdm, 0, &comm_1d_x, 111, 222);

    if (!is_serial[1])
        geometry_exchange_one_direction_gpu(
            u, sdm, 1, &comm_1d_y, 333, 444);

    if (!is_serial[2])
        geometry_exchange_one_direction_gpu(
            u, sdm, 2, &comm_1d_z, 555, 666);
}
