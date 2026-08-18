#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>

#define CUDA_CHECK(call) do {                                      \
    cudaError_t err = call;                                        \
    if (err != cudaSuccess) {                                      \
        printf("CUDA error: %s\n", cudaGetErrorString(err));       \
        exit(1);                                                   \
    }                                                             \
} while (0)

#ifndef COEFF_GPU
#define COEFF_GPU(a, m, i, j, k) \
    ((a)->coeff + (((m) * ((a)->nx + 1) * ((a)->ny + 1) * ((a)->nz + 1)) \
                 + ((i) * ((a)->ny + 1) * ((a)->nz + 1)) \
                 + ((j) * ((a)->nz + 1)) \
                 + (k)))
#endif

__global__
void matrix_poisson_create_kernel(matrix_poisson a_poisson, const double *dxm,const double *dym,const double *dzm,const double *dxg,const double *dyg,const double *dzg)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

    if (i > a_poisson.nx || j > a_poisson.ny || k > a_poisson.nz) return;

    double dxmp2i = 1.0 / (dxm[i] * dxg[i]);
    double dxmn2i = 1.0 / (dxm[i] * dxg[i+1]);

    double dymp2i = 1.0 / (dym[j] * dyg[j]);
    double dymn2i = 1.0 / (dym[j] * dyg[j+1]);

    double dzmp2i = 1.0 / (dzm[k] * dzg[k]);
    double dzmn2i = 1.0 / (dzm[k] * dzg[k+1]);

    *COEFF_GPU(&a_poisson, 0, i, j, k) = -(dxmp2i + dxmn2i + dymp2i + dymn2i + dzmp2i + dzmn2i);

    *COEFF_GPU(&a_poisson, 1, i, j, k) = dxmp2i;
    *COEFF_GPU(&a_poisson, 2, i, j, k) = dxmn2i;
    *COEFF_GPU(&a_poisson, 3, i, j, k) = dymp2i;
    *COEFF_GPU(&a_poisson, 4, i, j, k) = dymn2i;
    *COEFF_GPU(&a_poisson, 5, i, j, k) = dzmp2i;
    *COEFF_GPU(&a_poisson, 6, i, j, k) = dzmn2i;
}

void matrix_poisson_create_gpu(matrix_poisson *a_poisson, const subdomain *sdm)
{
    a_poisson->nx  = sdm->nx;
    a_poisson->ny  = sdm->ny;
    a_poisson->nz  = sdm->nz;
    a_poisson->dof = sdm->nx * sdm->ny * sdm->nz;

    size_t total_size = 7 * (a_poisson->nx + 1)
                          * (a_poisson->ny + 1)
                          * (a_poisson->nz + 1);

    CUDA_CHECK(cudaMalloc((void**)&a_poisson->coeff, total_size * sizeof(double)));

    CUDA_CHECK(cudaMemset(a_poisson->coeff, 0, total_size * sizeof(double)));

    dim3 block(8, 8, 4);

    dim3 grid((a_poisson->nx + block.x - 1) / block.x,
              (a_poisson->ny + block.y - 1) / block.y,
              (a_poisson->nz + block.z - 1) / block.z);

    matrix_poisson_create_kernel<<<grid, block>>>(*a_poisson, sdm->dxm, sdm->dym, sdm->dzm, sdm->dxg, sdm->dyg, sdm->dzg);

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
}

void matrix_poisson_destroy_gpu(matrix_poisson *a_poisson)
{
    a_poisson->dof = 0;

    if (a_poisson->coeff != NULL) {
        CUDA_CHECK(cudaFree(a_poisson->coeff));
        a_poisson->coeff = NULL;
    }
}