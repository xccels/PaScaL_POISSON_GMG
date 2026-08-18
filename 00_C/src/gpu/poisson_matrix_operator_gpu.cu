#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "poisson_matrix_operator.h"
// #include "timer.h"

#define CUDA_CHECK(call)                                      \
do {                                                          \
    cudaError_t err = call;                                   \
    if (err != cudaSuccess) {                                 \
        printf("CUDA error at %s:%d: %s\n",                  \
               __FILE__, __LINE__, cudaGetErrorString(err)); \
        exit(1);                                              \
    }                                                         \
} while (0)

/* ============================================================
   Dot product:
   result = sum x[idx] * y[idx], interior cells only
   ============================================================ */

__global__
void vv_dot_3d_matrix_kernel(const double *x, const double *y, double *partial_sum, int nx, int ny, int nz)
{
    extern __shared__ double sdata[];

    int tid = threadIdx.x;
    int gid = blockIdx.x * blockDim.x + threadIdx.x;

    int ni = (ny + 2) * (nz + 2);
    int nj = nz + 2;

    int total_inner = nx * ny * nz;

    double sum = 0.0;

    for (int id = gid; id < total_inner; id += blockDim.x * gridDim.x)
    {
        int kk = id % nz;
        int jj = (id / nz) % ny;
        int ii = id / (ny * nz);

        int i = ii + 1;
        int j = jj + 1;
        int k = kk + 1;

        size_t idx = (size_t)i * ni + (size_t)j * nj + k;

        sum += x[idx] * y[idx];
    }

    sdata[tid] = sum;
    __syncthreads();

    for (int stride = blockDim.x / 2; stride > 0; stride >>= 1)
    {
        if (tid < stride)
        {
            sdata[tid] += sdata[tid + stride];
        }
        __syncthreads();
    }

    if (tid == 0)
    {
        partial_sum[blockIdx.x] = sdata[0];
    }
}

void vv_dot_3d_matrix_gpu(double *result, const double *d_x, const double *d_y, int nx, int ny, int nz, int is_serial[3])
{
    int threads = 256;
    int total_inner = nx * ny * nz;

    int blocks = (total_inner + threads - 1) / threads;
    if (blocks < 1) blocks = 1;
    if (blocks > 256) blocks = 256;

    double *d_partial_sum = NULL;
    double *h_partial_sum = NULL;

    h_partial_sum = (double*)malloc(blocks * sizeof(double));
    if (h_partial_sum == NULL) {
        printf("Host malloc failed in vv_dot_3d_matrix_gpu\n");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    CUDA_CHECK(cudaMalloc((void**)&d_partial_sum, blocks * sizeof(double)));

    // timer_stamp0(STAMP_COMP);
    vv_dot_3d_matrix_kernel<<<blocks, threads, threads * sizeof(double)>>>(d_x, d_y, d_partial_sum, nx, ny, nz);
    
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    // timer_stamp(10,STAMP_COMP);

    CUDA_CHECK(cudaMemcpy(h_partial_sum, d_partial_sum, blocks * sizeof(double), cudaMemcpyDeviceToHost));

    // timer_stamp0(STAMP_COMP);
    double result_local = 0.0;
    for (int i = 0; i < blocks; i++) {
        result_local += h_partial_sum[i];
    }
    // timer_stamp(10,STAMP_COMP);

    double result_x, result_xy;

    if (is_serial[0] && is_serial[1] && is_serial[2]) {
        // timer_stamp0(STAMP_COMP);
        *result = result_local;
        // timer_stamp(10,STAMP_COMP);

    } else if (!is_serial[0] && !is_serial[1] && !is_serial[2]) {
        // timer_stamp0(STAMP_COMM_ALLREDUCE);
        MPI_Allreduce(&result_local, result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        // timer_stamp(12,STAMP_COMM_ALLREDUCE);

    } else {
        if (!is_serial[0]) {
            // timer_stamp0(STAMP_COMM_ALLREDUCE);
            MPI_Allreduce(&result_local, &result_x, 1, MPI_DOUBLE, MPI_SUM, comm_1d_x.mpi_comm);
            // timer_stamp(12,STAMP_COMM_ALLREDUCE);
        } else {
            // timer_stamp0(STAMP_COMP);
            result_x = result_local;
            // timer_stamp(10,STAMP_COMP);
        }

        if (!is_serial[1]) {
            // timer_stamp0(STAMP_COMM_ALLREDUCE);
            MPI_Allreduce(&result_x, &result_xy, 1, MPI_DOUBLE, MPI_SUM, comm_1d_y.mpi_comm);
            // timer_stamp(12,STAMP_COMM_ALLREDUCE);
        } else {
            // timer_stamp0(STAMP_COMP);
            result_xy = result_x;
            // timer_stamp(10,STAMP_COMP);
        }

        if (!is_serial[2]) {
            // timer_stamp0(STAMP_COMM_ALLREDUCE);
            MPI_Allreduce(&result_xy, result, 1, MPI_DOUBLE, MPI_SUM, comm_1d_z.mpi_comm);
            // timer_stamp(12,STAMP_COMM_ALLREDUCE);

        } else {
            // timer_stamp0(STAMP_COMP);
            *result = result_xy;
            // timer_stamp(10,STAMP_COMP);
        }
    }

    CUDA_CHECK(cudaFree(d_partial_sum));
    free(h_partial_sum);
}

/* ============================================================
   Set full 3D array to zero, including halo cells
   ============================================================ */

__global__
void set_zero_3d_kernel(double *y, int nx, int ny, int nz)
{
    int gid = blockIdx.x * blockDim.x + threadIdx.x;

    int total_all = (nx + 2) * (ny + 2) * (nz + 2);

    for (int id = gid; id < total_all; id += blockDim.x * gridDim.x)
    {
        y[id] = 0.0;
    }
}

/* ============================================================
   Poisson matrix-vector multiplication:
   y = A x

   x/y layout: (nx+2) * (ny+2) * (nz+2)

   coef layout: 7 * (nx+1) * (ny+1) * (nz+1)

   coef[0] : center
   coef[1] : i-1
   coef[2] : i+1
   coef[3] : j-1
   coef[4] : j+1
   coef[5] : k-1
   coef[6] : k+1
   ============================================================ */

__global__
void mv_mul_poisson_matrix_kernel(double *y, const double *x, const double *coef, int nx, int ny, int nz)
{
    int gid = blockIdx.x * blockDim.x + threadIdx.x;

    int ni  = (ny + 2) * (nz + 2);
    int nj  = nz + 2;

    int cni = (ny + 1) * (nz + 1);
    int cnj = nz + 1;

    int coef_size = (nx + 1) * (ny + 1) * (nz + 1);

    int total_inner = nx * ny * nz;

    for (int id = gid; id < total_inner; id += blockDim.x * gridDim.x)
    {
        int kk = id % nz;
        int jj = (id / nz) % ny;
        int ii = id / (ny * nz);

        int i = ii + 1;
        int j = jj + 1;
        int k = kk + 1;

        size_t idx  = (size_t)i * ni  + (size_t)j * nj  + k;
        size_t cidx = (size_t)i * cni + (size_t)j * cnj + k;

        y[idx] =
              coef[0 * coef_size + cidx] * x[idx]
            + coef[1 * coef_size + cidx] * x[idx - ni]
            + coef[2 * coef_size + cidx] * x[idx + ni]
            + coef[3 * coef_size + cidx] * x[idx - nj]
            + coef[4 * coef_size + cidx] * x[idx + nj]
            + coef[5 * coef_size + cidx] * x[idx - 1]
            + coef[6 * coef_size + cidx] * x[idx + 1];
    }
}

void mv_mul_poisson_matrix_gpu(double *d_y, const double *d_coef, double *d_x, subdomain *dm, int is_serial[3])
{
    int nx = dm->nx;
    int ny = dm->ny;
    int nz = dm->nz;

    if (!(is_serial[0] && is_serial[1] && is_serial[2])) {
        CUDA_CHECK(cudaDeviceSynchronize());

        // timer_stamp0(STAMP_COMM_NEIGHBOR);
        geometry_halocell_update_selectively_gpu(d_x, dm, is_serial);
        // timer_stamp(11,STAMP_COMM_NEIGHBOR);

        CUDA_CHECK(cudaDeviceSynchronize());
    }

    int threads = 256;

    int total_all   = (nx + 2) * (ny + 2) * (nz + 2);
    int total_inner = nx * ny * nz;

    int blocks_zero  = (total_all + threads - 1) / threads;
    int blocks_inner = (total_inner + threads - 1) / threads;

    if (blocks_zero < 1)  blocks_zero = 1;
    if (blocks_inner < 1) blocks_inner = 1;

    // timer_stamp0(STAMP_COMP);
    set_zero_3d_kernel<<<blocks_zero, threads>>>(d_y, nx, ny, nz);
    // timer_stamp(10,STAMP_COMP);

    CUDA_CHECK(cudaGetLastError());

    // timer_stamp0(STAMP_COMP);
    mv_mul_poisson_matrix_kernel<<<blocks_inner, threads>>>(d_y, d_x, d_coef, nx, ny, nz);
    // timer_stamp(10,STAMP_COMP);
    
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
}
