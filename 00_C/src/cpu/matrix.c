#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>


void matrix_poisson_create(matrix_poisson *a_poisson, const subdomain *sdm) {
    int i, j, k;
    double dxm2i, dym2i, dzm2i;
    double dxmp2i, dxmn2i, dymp2i, dymn2i, dzmp2i, dzmn2i;
    
    a_poisson->nx = sdm->nx;
    a_poisson->ny = sdm->ny;
    a_poisson->nz = sdm->nz;
    a_poisson->dof = sdm->nx * sdm->ny * sdm->nz;

    size_t total_size = 7 * (a_poisson->nx + 1) * (a_poisson->ny + 1) * (a_poisson->nz + 1);
    a_poisson->coeff = calloc(total_size, sizeof(double));
    if (!a_poisson->coeff) {
        fprintf(stderr, "Error: memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    for (i = 1; i <= sdm->nx; i++) {
        dxm2i = 1.0 / (sdm->dxm[i] * sdm->dxm[i]);
        dxmp2i = 1.0 / (sdm->dxm[i] * sdm->dxg[i]);
        dxmn2i = 1.0 / (sdm->dxm[i] * sdm->dxg[i+1]);

        for (j = 1; j <= sdm->ny; j++) {
            dym2i = 1.0 / (sdm->dym[j] * sdm->dym[j]);
            dymp2i = 1.0 / (sdm->dym[j] * sdm->dyg[j]);
            dymn2i = 1.0 / (sdm->dym[j] * sdm->dyg[j+1]);

            for (k = 1; k <= sdm->nz; k++) {
                dzm2i = 1.0 / (sdm->dzm[k] * sdm->dzm[k]);
                dzmp2i = 1.0 / (sdm->dzm[k] * sdm->dzg[k]);
                dzmn2i = 1.0 / (sdm->dzm[k] * sdm->dzg[k+1]);

                *COEFF(a_poisson, 0, i, j, k) = -(dxmp2i + dxmn2i + dymp2i + dymn2i + dzmp2i + dzmn2i);
                *COEFF(a_poisson, 1, i, j, k) = dxmp2i;
                *COEFF(a_poisson, 2, i, j, k) = dxmn2i;
                *COEFF(a_poisson, 3, i, j, k) = dymp2i;
                *COEFF(a_poisson, 4, i, j, k) = dymn2i;
                *COEFF(a_poisson, 5, i, j, k) = dzmp2i;
                *COEFF(a_poisson, 6, i, j, k) = dzmn2i;
            }
        }
    }
}

void matrix_poisson_destroy(matrix_poisson *a_poisson) 
{
    a_poisson->dof = 0;

    if (a_poisson->coeff != NULL) {
        free(a_poisson->coeff);
        a_poisson->coeff = NULL;
    }
}