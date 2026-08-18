// poisson_matrix_operator.c
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdbool.h>
#include <stddef.h>
#include "poisson_matrix_operator.h"
#include "timer.h"

void vv_dot_3d_matrix(double *result, double *x, double *y, int nx, int ny, int nz, int is_serial[3]) 
{
    int i, j, k, ni,nj;
    double result_local = 0.0;
    double result_x, result_xy;

    timer_stamp0(STAMP_COMP);
    ni = (ny+2)*(nz+2);
    nj = (nz+2);
    #pragma omp parallel for collapse(3) reduction(+:result_local)
    for(i=1; i<=nx; i++){
        for(j=1; j<=ny; j++){
            for(k=1; k<=nz; k++){
                size_t idx = (i)*ni + (j)*nj + (k);
                result_local += x[idx] * y[idx];
            }
        }
    }
    timer_stamp(10,STAMP_COMP);


    if (is_serial[0] && is_serial[1] && is_serial[2]) {
        timer_stamp0(STAMP_COMP);
        *result = result_local;
        timer_stamp(10,STAMP_COMP);

    } else if (!is_serial[0] && !is_serial[1] && !is_serial[2]) {
        timer_stamp0(STAMP_COMM_ALLREDUCE);
        MPI_Allreduce(&result_local, result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        timer_stamp(12,STAMP_COMM_ALLREDUCE);

    } else {
        if (!is_serial[0]) {
            timer_stamp0(STAMP_COMM_ALLREDUCE);
            MPI_Allreduce(&result_local, &result_x, 1, MPI_DOUBLE, MPI_SUM, comm_1d_x.mpi_comm);
            timer_stamp(12,STAMP_COMM_ALLREDUCE);

        } else {
            timer_stamp0(STAMP_COMP);
            result_x = result_local;
            timer_stamp(10,STAMP_COMP);

        }
        if (!is_serial[1]) {
            timer_stamp0(STAMP_COMM_ALLREDUCE);
            MPI_Allreduce(&result_x, &result_xy, 1, MPI_DOUBLE, MPI_SUM, comm_1d_y.mpi_comm);
            timer_stamp(12,STAMP_COMM_ALLREDUCE);

        } else {
            timer_stamp0(STAMP_COMP);
            result_xy = result_x;
            timer_stamp(10,STAMP_COMP);

        }
        if (!is_serial[2]) {
            timer_stamp0(STAMP_COMM_ALLREDUCE);
            MPI_Allreduce(&result_xy, result, 1, MPI_DOUBLE, MPI_SUM, comm_1d_z.mpi_comm);
            timer_stamp(12,STAMP_COMM_ALLREDUCE);

        } else {
            timer_stamp0(STAMP_COMP);
            *result = result_xy;
            timer_stamp(10,STAMP_COMP);
        }
    }
}

void mv_mul_poisson_matrix(double *y, matrix_poisson *a_poisson, double *x, subdomain *dm, int is_serial[3]) 
{
    
    int i, j, k, ni, nj;
    
    int nx = dm->nx;
    int ny = dm->ny;
    int nz = dm->nz;

    // printf("[Poisson_matrix_operator] is_serial[0] = %d\n", is_serial[0]);
    // printf("[Poisson_matrix_operator] is_serial[1] = %d\n", is_serial[1]);
    // printf("[Poisson_matrix_operator] is_serial[2] = %d\n", is_serial[2]);
    if (!(is_serial[0] && is_serial[1] && is_serial[2])) {
        timer_stamp0(STAMP_COMM_NEIGHBOR);
        geometry_halocell_update_selectively(x, dm, is_serial);
        timer_stamp(11,STAMP_COMM_NEIGHBOR);

    }

    timer_stamp0(STAMP_COMP);
    ni = (ny+2)*(nz+2);
    nj = (nz+2);
    #pragma omp parallel for
    for(i=0; i<=nx+1; i++){
        for(j=0; j<=ny+1; j++){
            for(k=0; k<=nz+1; k++){
                size_t idx = (i)*ni + (j)*nj + (k);
                y[idx] = 0.0;
            }
        }
    }

    #pragma omp parallel for
    for(i=1; i<=nx; i++){
        for(j=1; j<=ny; j++){
            for(k=1; k<=nz; k++){
                size_t idx = (i)*ni + (j)*nj + (k);
                y[idx] = *COEFF(a_poisson, 0, i, j, k) * x[idx] 
                       + *COEFF(a_poisson, 1, i, j, k) * x[idx - ni] 
                       + *COEFF(a_poisson, 2, i, j, k) * x[idx + ni] 
                       + *COEFF(a_poisson, 3, i, j, k) * x[idx - nj] 
                       + *COEFF(a_poisson, 4, i, j, k) * x[idx + nj] 
                       + *COEFF(a_poisson, 5, i, j, k) * x[idx - 1] 
                       + *COEFF(a_poisson, 6, i, j, k) * x[idx + 1];
            }
        }
    }
    timer_stamp(10,STAMP_COMP);

    // timer_stamp(10, 10);
}
