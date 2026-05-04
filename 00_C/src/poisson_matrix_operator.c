// poisson_matrix_operator.c
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdbool.h>
#include <stddef.h>
#include "poisson_matrix_operator.h"

void vv_dot_3d_matrix(double *result, double *x, double *y, int nx, int ny, int nz, int is_serial[3]) 
{
    int i, j, k;
    double result_local = 0.0;
    double result_x, result_xy;

    // timer_stamp0(10); // stamp_comp

    #pragma omp parallel for collapse(3) reduction(+:result_local)
    for(i=1; i<=nx; i++){
        for(j=1; j<=ny; j++){
            for(k=1; k<=nz; k++){
                size_t idx = (i)*(ny+2)*(nz+2) + (j)*(nz+2) + (k);
                result_local += x[idx] * y[idx];
            }
        }
    }

    // timer_stamp(10, 10); // stamp_comp

    if (is_serial[0] && is_serial[1] && is_serial[2]) {
        // timer_stamp0(10);
        *result = result_local;
        // timer_stamp(10, 10);
    } else if (!is_serial[0] && !is_serial[1] && !is_serial[2]) {
        // timer_stamp0(12);
        MPI_Allreduce(&result_local, result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        // timer_stamp(12, 12);
    } else {
        if (!is_serial[0]) {
            // timer_stamp0(12);
            MPI_Allreduce(&result_local, &result_x, 1, MPI_DOUBLE, MPI_SUM, comm_1d_x.mpi_comm);
            // timer_stamp(12, 12);
        } else {
            // timer_stamp0(10);
            result_x = result_local;
            // timer_stamp(10, 10);
        }
        if (!is_serial[1]) {
            // timer_stamp0(12);
            MPI_Allreduce(&result_x, &result_xy, 1, MPI_DOUBLE, MPI_SUM, comm_1d_y.mpi_comm);
            // timer_stamp(12, 12);
        } else {
            // timer_stamp0(10);
            result_xy = result_x;
            // timer_stamp(10, 10);
        }
        if (!is_serial[2]) {
            // timer_stamp0(12);
            MPI_Allreduce(&result_xy, result, 1, MPI_DOUBLE, MPI_SUM, comm_1d_z.mpi_comm);
            // timer_stamp(12, 12);
        } else {
            // timer_stamp0(10);
            *result = result_xy;
            // timer_stamp(10, 10);
        }
    }
}

void mv_mul_poisson_matrix(double *y, matrix_poisson *a_poisson, double *x, subdomain *dm, int is_serial[3]) 
{
    
    int i, j, k;
    
    int nx = dm->nx;
    int ny = dm->ny;
    int nz = dm->nz;

    // printf("[Poisson_matrix_operator] is_serial[0] = %d\n", is_serial[0]);
    // printf("[Poisson_matrix_operator] is_serial[1] = %d\n", is_serial[1]);
    // printf("[Poisson_matrix_operator] is_serial[2] = %d\n", is_serial[2]);
    if (!(is_serial[0] && is_serial[1] && is_serial[2])) {
        // timer_stamp0(11);
        geometry_halocell_update_selectively(x, dm, is_serial);
        // printf("[Poisson_matrix_operator] here\n");
        // timer_stamp(11, 11);
    }

    // timer_stamp0(10);

    #pragma omp parallel for
    for(i=0; i<=nx+1; i++){
        for(j=0; j<=ny+1; j++){
            for(k=0; k<=nz+1; k++){
                size_t idx = (i)*(ny+2)*(nz+2) + (j)*(nz+2) + (k);
                y[idx] = 0.0;
            }
        }
    }

    #pragma omp parallel for
    for(i=1; i<=nx; i++){
        for(j=1; j<=ny; j++){
            for(k=1; k<=nz; k++){
                size_t idx = (i)*(ny+2)*(nz+2) + (j)*(nz+2) + (k);
                y[idx] = *COEFF(a_poisson, 0, i, j, k) * x[idx] 
                       + *COEFF(a_poisson, 1, i, j, k) * x[idx - (ny+2)*(nz+2)] 
                       + *COEFF(a_poisson, 2, i, j, k) * x[idx + (ny+2)*(nz+2)] 
                       + *COEFF(a_poisson, 3, i, j, k) * x[idx - (nz+2)] 
                       + *COEFF(a_poisson, 4, i, j, k) * x[idx + (nz+2)] 
                       + *COEFF(a_poisson, 5, i, j, k) * x[idx - 1] 
                       + *COEFF(a_poisson, 6, i, j, k) * x[idx + 1];
            }
        }
    }

    // timer_stamp(10, 10);
}
