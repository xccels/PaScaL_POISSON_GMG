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
#include "global.h"


void rbgs_solver_poisson_matrix(double *sol, matrix_poisson *a_poisson, double *rhs, subdomain *dm, int maxiteration, double tolerance, double omega, int is_aggregated[3]);

void rbgs_iterator_poisson_matrix(double *sol, matrix_poisson *a_poisson, double *rhs, subdomain *dm, int maxiteration, double omega, int is_aggregated[3]);

#endif // RBGS_POISSON_MATRIX_H
