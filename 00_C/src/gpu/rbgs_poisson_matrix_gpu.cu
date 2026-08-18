#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda_runtime.h>

#include "rbgs_poisson_matrix.h"
#include "geometry.h"
// #include "timer.h" 

#ifndef IDX
#define IDX(i,j,k,ni,nj) ((i)*(ni) + (j)*(nj) + (k))
#endif

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
   rsd = rhs - rsd
   ============================================================ */

__global__
static void compute_residual_kernel(double *rsd, const double *rhs, int nx, int ny, int nz)
{
    int gid = blockIdx.x * blockDim.x + threadIdx.x;

    int ni = (ny + 2) * (nz + 2);
    int nj = nz + 2;

    int total_inner = nx * ny * nz;

    for (int id = gid; id < total_inner; id += blockDim.x * gridDim.x)
    {
        int kk = id % nz;
        int jj = (id / nz) % ny;
        int ii = id / (ny * nz);

        int i = ii + 1;
        int j = jj + 1;
        int k = kk + 1;

        size_t idx = (size_t)i * ni + (size_t)j * nj + k;

        rsd[idx] = rhs[idx] - rsd[idx];
    }
}

/* ============================================================
   residual norm
   rsd = rhs - A*x
   ============================================================ */

void compute_residual_norm_gpu(double *rsd_norm, double *d_rsd, double *d_coef, double *d_sol, double *d_rhs, subdomain *dm, int is_aggregated[3])
{
    int nx = dm->nx;
    int ny = dm->ny;
    int nz = dm->nz;

    int threads = 256;
    int total_inner = nx * ny * nz;
    int blocks = (total_inner + threads - 1) / threads;

    CUDA_CHECK(cudaMemset(d_rsd, 0, (size_t)(nx + 2) * (ny + 2) * (nz + 2) * sizeof(double)));

    mv_mul_poisson_matrix_gpu(d_rsd, d_coef, d_sol, dm, is_aggregated);

    // timer_stamp0(STAMP_COMP);
    compute_residual_kernel<<<blocks, threads>>>(d_rsd, d_rhs, nx, ny, nz);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    // timer_stamp(10,STAMP_COMP);

    vv_dot_3d_matrix_gpu(rsd_norm, d_rsd, d_rsd, nx, ny, nz, is_aggregated);
}

/* ============================================================
   One red-black Gauss-Seidel sweep
   color = 0 or 1

   sol/rhs index: (nx+2)(ny+2)(nz+2)

   coef index: (nx+1)(ny+1)(nz+1)
   ============================================================ */

__global__
static void rbgs_sweep_kernel(double *sol, const double *rhs, const double *coef, int nx, int ny, int nz, double omega, int sweep, int offset)
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

        int k_start = 1 + ((i + j + sweep + offset) % 2);
        if (((k - k_start) % 2) != 0) continue;

        size_t idx  = (size_t)i * ni  + (size_t)j * nj  + k;
        size_t cidx = (size_t)i * cni + (size_t)j * cnj + k;

        double ac = coef[0 * coef_size + cidx];

        double temp =
              coef[1 * coef_size + cidx] * sol[idx - ni]
            + coef[2 * coef_size + cidx] * sol[idx + ni]
            + coef[3 * coef_size + cidx] * sol[idx - nj]
            + coef[4 * coef_size + cidx] * sol[idx + nj]
            + coef[5 * coef_size + cidx] * sol[idx - 1]
            + coef[6 * coef_size + cidx] * sol[idx + 1];

        sol[idx] = omega * (rhs[idx] - temp) / ac
                 + (1.0 - omega) * sol[idx];
    }
}

/* ============================================================
   GPU RBGS solver
   with tolerance check
   ============================================================ */

void rbgs_solver_poisson_matrix_gpu(double *d_sol, double *d_coef, double *d_rhs, subdomain *dm, int maxiteration, double tolerance, double omega, int is_aggregated[3])
{
    int nx = dm->nx;
    int ny = dm->ny;
    int nz = dm->nz;

    size_t total_size =
        (size_t)(nx + 2) *
        (size_t)(ny + 2) *
        (size_t)(nz + 2);

    int total_inner = nx * ny * nz;

    int threads = 256;
    int blocks = (total_inner + threads - 1) / threads;

    double *d_rsd = NULL;
    CUDA_CHECK(cudaMalloc((void**)&d_rsd, total_size * sizeof(double)));

    double rsd0tol = 0.0;
    double rsd_norm = 0.0;

    compute_residual_norm_gpu(&rsd0tol, d_rsd, d_coef, d_sol, d_rhs, dm, is_aggregated);

    // timer_stamp0(STAMP_COMP);
    rsd_norm = rsd0tol;

    int is_same[3];

    is_same[0] = is_aggregated[0] ? 0 : 1;
    is_same[1] = is_aggregated[1] ? 0 : 1;
    is_same[2] = is_aggregated[2] ? 0 : 1;

    int offset = ( (nx % 2) * ((comm_1d_x.myrank * is_same[0]) % 2)
                 + (ny % 2) * ((comm_1d_y.myrank * is_same[1]) % 2)
                 + (nz % 2) * ((comm_1d_z.myrank * is_same[2]) % 2) ) % 2;
    // timer_stamp(10,STAMP_COMP);

    int iter;

    for (iter = 0; iter < maxiteration; iter++)
    {
        // if ((iter % 10 == 0) && (myrank == 0)) {
        //     printf("[RBGS GPU] iter = %d, mse = %e, r_mse = %e\n", iter, sqrt(rsd_norm), sqrt(rsd_norm / rsd0tol));
        // }


        if (!(is_aggregated[0] && is_aggregated[1] && is_aggregated[2])) {
            // timer_stamp0(STAMP_COMM_NEIGHBOR);
            geometry_halocell_update_selectively_gpu(d_sol, dm, is_aggregated);
            // timer_stamp(11,STAMP_COMM_NEIGHBOR);
        }

        // timer_stamp0(STAMP_COMP);
        rbgs_sweep_kernel<<<blocks, threads>>>(d_sol, d_rhs, d_coef, nx, ny, nz, omega, 0, offset);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
        // timer_stamp(10,STAMP_COMP);

        if (!(is_aggregated[0] && is_aggregated[1] && is_aggregated[2])) {
            // timer_stamp0(STAMP_COMM_NEIGHBOR);
            geometry_halocell_update_selectively_gpu(d_sol, dm, is_aggregated);
            // timer_stamp(11,STAMP_COMM_NEIGHBOR);
        }

        // timer_stamp0(STAMP_COMP);
        rbgs_sweep_kernel<<<blocks, threads>>>(d_sol, d_rhs, d_coef, nx, ny, nz, omega, 1, offset); 
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
        // timer_stamp(10,STAMP_COMP);


        compute_residual_norm_gpu(&rsd_norm, d_rsd, d_coef, d_sol, d_rhs, dm, is_aggregated);

        if (sqrt(rsd_norm / rsd0tol) <= tolerance) {
            break;
        }
    }

    // if (myrank == 0) {
    //     printf("[RBGS GPU] completed. Iteration = %d, mse = %e, r_mse = %e\n", iter, sqrt(rsd_norm), sqrt(rsd_norm / rsd0tol));
    // }
        
    CUDA_CHECK(cudaFree(d_rsd));
}

/* ============================================================
   GPU RBGS iterator
   fixed iteration count
   ============================================================ */

void rbgs_iterator_poisson_matrix_gpu(double *d_sol, double *d_coef, double *d_rhs, subdomain *dm, int maxiteration, double omega, int is_aggregated[3])
{
    int nx = dm->nx;
    int ny = dm->ny;
    int nz = dm->nz;

    // timer_stamp0(STAMP_COMP);
    size_t total_size =
        (size_t)(nx + 2) *
        (size_t)(ny + 2) *
        (size_t)(nz + 2);

    int total_inner = nx * ny * nz;

    int threads = 256;
    int blocks = (total_inner + threads - 1) / threads;
    // timer_stamp(10,STAMP_COMP);

    double *d_rsd = NULL;
    CUDA_CHECK(cudaMalloc((void**)&d_rsd, total_size * sizeof(double)));

    // timer_stamp0(STAMP_COMP);
    double rsd0tol = 0.0;
    double rsd_norm = 0.0;

    int is_same[3];

    is_same[0] = is_aggregated[0] ? 0 : 1;
    is_same[1] = is_aggregated[1] ? 0 : 1;
    is_same[2] = is_aggregated[2] ? 0 : 1;

    int offset = ( (nx % 2) * ((comm_1d_x.myrank * is_same[0]) % 2)
                 + (ny % 2) * ((comm_1d_y.myrank * is_same[1]) % 2)
                 + (nz % 2) * ((comm_1d_z.myrank * is_same[2]) % 2) ) % 2;
    // timer_stamp(10,STAMP_COMP);

    compute_residual_norm_gpu(&rsd0tol, d_rsd, d_coef, d_sol, d_rhs, dm, is_aggregated);

    for (int iter = 0; iter < maxiteration; iter++)
    {
        if (!(is_aggregated[0] && is_aggregated[1] && is_aggregated[2])) {
            // timer_stamp0(STAMP_COMM_NEIGHBOR);
            geometry_halocell_update_selectively_gpu(d_sol, dm, is_aggregated);
            // timer_stamp(11,STAMP_COMM_NEIGHBOR);
        }

        // timer_stamp0(STAMP_COMP);
        rbgs_sweep_kernel<<<blocks, threads>>>(d_sol, d_rhs, d_coef, nx, ny, nz, omega, 0, offset);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
        // timer_stamp(10,STAMP_COMP);

        if (!(is_aggregated[0] && is_aggregated[1] && is_aggregated[2])) {
            // timer_stamp0(STAMP_COMM_NEIGHBOR);
            geometry_halocell_update_selectively_gpu(d_sol, dm, is_aggregated);
            // timer_stamp(11,STAMP_COMM_NEIGHBOR);
        }

        // timer_stamp0(STAMP_COMP);
        rbgs_sweep_kernel<<<blocks, threads>>>(d_sol, d_rhs, d_coef, nx, ny, nz, omega, 1, offset);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
        // timer_stamp(10,STAMP_COMP);

    }

    compute_residual_norm_gpu(&rsd_norm, d_rsd, d_coef, d_sol, d_rhs, dm, is_aggregated);

    // if (myrank == 0) {
    //     printf("[RBGS iterator GPU] completed. Iteration = %d, mse = %e, r_mse = %e, rsd0 = %e \n", maxiteration, sqrt(rsd_norm), sqrt(rsd_norm/rsd0tol), sqrt(rsd0tol));
    // }
    
    CUDA_CHECK(cudaFree(d_rsd));
}



