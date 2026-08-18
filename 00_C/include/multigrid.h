#ifndef MULTIGRID_H
#define MULTIGRID_H

#include "geometry.h"  // subdomain 类型
#include "matrix.h"    // matrix_poisson 类型
#include "mpi_topology.h"   
#include "rbgs_poisson_matrix.h" 
// #include "global.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================
   CPU interfaces
   ============================================================ */
void multigrid_create(subdomain *sdm, int nlevel, int ncycle, int aggr_method, int aggr_level);
void multigrid_solve_vcycle(double *sol, matrix_poisson *a_poisson, double *rhs, subdomain *sdm, int maxiteration, double tolerance, double omega_sor);
void multigrid_destroy(void);

void multigrid_subdomain_create_no_aggregation(subdomain *sdm);
void multigrid_subdomain_create_single_aggregation(subdomain *sdm);
void multigrid_subdomain_create_adaptive_aggregation(subdomain *sdm);
void multigrid_subdomain_no_aggregation_make_grid(double *dxm, double *dxg, double *xg, int nx,
                                                  const double *dxm_f, const double *dxg_f, const double *xg_f,
                                                  double ox, int lv_cur, int lv_coarsest,
                                                  cart_comm_1d comm_1d, char dir);
void multigrid_subdomain_aggregation_make_grid(double *dxm, double *dxg, double *xg, int nx,
                                               const double *dxm_f, const double *dxg_f, const double *xg_f,
                                               double ox, int lv_cur, int lv_coarsest, int lv_aggregation, int nx_lv_aggregation,
                                               cart_comm_1d comm_1d, char dir);
void multigrid_allocate_subdomain_variables(void);

void multigrid_no_aggregation_vcycle_solver(double *sol, matrix_poisson *a_poisson, double *rhs, subdomain *sdm, int maxiteration, double tolerance, double omega_sor);
void multigrid_single_aggregation_vcycle_solver(double *sol, matrix_poisson *a_poisson, double *rhs, subdomain *sdm, int maxiteration, double tolerance, double omega_sor);
void multigrid_adaptive_aggregation_vcycle_solver(double *sol, matrix_poisson *a_poisson, double *rhs, subdomain *sdm, int maxiteration, double tolerance, double omega_sor);
void multigrid_solve_coarset_level(double *x, matrix_poisson *a_poisson, double *rhs, subdomain *dm, int maxiteration, double tolerance, double omega, int is_aggregated[3]);

/* ============================================================
   GPU interfaces
   ============================================================ */
   void multigrid_create_gpu(subdomain *sdm, int nlevel, int ncycle, int aggr_method, int aggr_level);
   void multigrid_subdomain_create_no_aggregation_gpu(subdomain *sdm);
   void multigrid_subdomain_create_single_aggregation_gpu(subdomain *sdm);
   void multigrid_subdomain_create_adaptive_aggregation_gpu(subdomain *sdm);
   void multigrid_subdomain_no_aggregation_make_grid_gpu( double *dxm, double *dxg, double *xg, int nx, const double *dxm_f, const double *dxg_f, const double *xg_f, double ox, int lv_cur, int lv_coarsest, cart_comm_1d comm_1d, char dir);
   void multigrid_subdomain_aggregation_make_grid_gpu(double *dxm, double *dxg, double *xg, int nx, const double *dxm_f, const double *dxg_f, const double *xg_f, double ox, int lv_cur, int lv_coarsest, int lv_aggregation, int nx_lv_aggregation, cart_comm_1d comm_1d, char dir);
   void multigrid_allocate_subdomain_variables_gpu(void);

   void multigrid_restriction_gpu(double *d_val_c, double *d_val_f, subdomain *dm_c, subdomain *dm_f, int level);
   void multigrid_prolongation_linear_on_nonuniform_grid_gpu(double *d_val_f, const double *d_val_c, subdomain *dm_f, subdomain *dm_c, int level);
   
   void multigrid_residual_gpu(double *d_rsd, double *d_coef, double *d_x, double *d_rhs, subdomain *dm, int is_aggregated[3]);
   void multigrid_solve_coarset_level_gpu(double *x, matrix_poisson *a_poisson, double *rhs, subdomain *dm, int maxiteration, double tolerance, double omega, int is_aggregated[3]);
   void multigrid_solve_vcycle_gpu(double *d_sol, matrix_poisson *a_poisson, double *d_rhs, subdomain *sdm, int maxiteration, double tolerance, double omega_sor);
   void multigrid_no_aggregation_vcycle_solver_gpu(double *sol, matrix_poisson *a_poisson, double *rhs, subdomain *sdm, int maxiteration, double tolerance, double omega_sor);
   void multigrid_single_aggregation_vcycle_solver_gpu(double *sol, matrix_poisson *a_poisson, double *rhs, subdomain *sdm, int maxiteration, double tolerance, double omega_sor);
   void multigrid_adaptive_aggregation_vcycle_solver_gpu(double *sol, matrix_poisson *a_poisson, double *rhs, subdomain *sdm, int maxiteration, double tolerance, double omega_sor);

   void multigrid_destroy_gpu();

#ifdef __cplusplus
}
#endif

#endif