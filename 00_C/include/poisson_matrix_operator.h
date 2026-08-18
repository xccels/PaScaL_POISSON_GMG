#ifndef POISSON_MATRIX_OPERATOR_H
#define POISSON_MATRIX_OPERATOR_H

#include "matrix.h"    // matrix_poisson
#include "geometry.h"  //subdomain
#include "timer.h"
#include "mpi_topology.h"
// #include "global.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================
   CPU interfaces
   ============================================================ */
   // 点积
   void vv_dot_3d_matrix(double *result, double *x, double *y, int nx, int ny, int nz, int is_serial[3]);

   // 矩阵乘法
   void mv_mul_poisson_matrix(double *y, matrix_poisson *a_poisson, double *x, subdomain *dm, int is_serial[3]);

/* ============================================================
   GPU interfaces
   ============================================================ */

   void vv_dot_3d_matrix_gpu(double *result, const double *d_x, const double *d_y, int nx, int ny, int nz, int is_serial[3]);

   void mv_mul_poisson_matrix_gpu(double *d_y, const double *d_coef, double *d_x, subdomain *dm, int is_serial[3]);
                               
#ifdef __cplusplus
}
#endif



#endif
