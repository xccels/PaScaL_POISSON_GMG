#ifndef MATRIX_H
#define MATRIX_H

#include "geometry.h"  // 引入 subdomain 结构体
#include <stdbool.h>
#include "global.h"

typedef struct matrix_poisson_struct {
    int nx, ny, nz;
    int dof;
    double *coeff;  // coeff[7*(nx+1)*(ny+1)*(nz+1)]
} matrix_poisson;

static inline double* COEFF(matrix_poisson *a, int m, int i, int j, int k) {
    return &a->coeff[ m*(a->nx+1)*(a->ny+1)*(a->nz+1)
                     + i*(a->ny+1)*(a->nz+1)
                     + j*(a->nz+1)
                     + k ];
}

void matrix_poisson_create(matrix_poisson *a_poisson, const subdomain *sdm);
void matrix_poisson_destroy(matrix_poisson *a_poisson);

#endif // MATRIX_H