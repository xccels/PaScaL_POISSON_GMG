#ifndef MULTIGRID_H
#define MULTIGRID_H

#include "geometry.h"  // subdomain 类型
#include "matrix.h"    // matrix_poisson 类型
#include "mpi_topology.h"   
#include "rbgs_poisson_matrix.h" 
#include "global.h"

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

#endif