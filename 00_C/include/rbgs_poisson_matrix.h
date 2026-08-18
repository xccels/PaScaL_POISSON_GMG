#ifndef RBGS_POISSON_MATRIX_H
#define RBGS_POISSON_MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>
#include "matrix.h"        // matrix_poisson
#include "geometry.h"      // subdomain, geometry_halocell_update_selectively
#include "mpi_topology.h"  // myrank, comm_1d_x, comm_1d_y, comm_1d_z
#include "poisson_matrix_operator.h" 
// #include "global.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================
   CPU interfaces
   ============================================================ */
void rbgs_solver_poisson_matrix(double *sol, matrix_poisson *a_poisson, double *rhs, subdomain *dm, int maxiteration, double tolerance, double omega, int is_aggregated[3]);

void rbgs_iterator_poisson_matrix(double *sol, matrix_poisson *a_poisson, double *rhs, subdomain *dm, int maxiteration, double omega, int is_aggregated[3]);


/* ============================================================
   GPU interfaces
   ============================================================ */

// void rbgs_gpu(double *d_sol, double *d_rhs, double *d_coef, int nx, int ny, int nz, int ni, int nj, double omega, int offset);
   void rbgs_solver_poisson_matrix_gpu(double *d_sol, double *d_coef, double *d_rhs, subdomain *dm, int maxiteration, double tolerance, double omega, int is_aggregated[3]);

   void rbgs_iterator_poisson_matrix_gpu(double *d_sol, double *d_coef, double *d_rhs, subdomain *dm, int maxiteration, double omega, int is_aggregated[3]);
   
#ifdef __cplusplus
}
#endif

#endif // RBGS_POISSON_MATRIX_H
