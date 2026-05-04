#include "rbgs_poisson_matrix.h" 

void rbgs_solver_poisson_matrix(double *sol, matrix_poisson *a_poisson, double *rhs, subdomain *dm, int maxiteration, double tolerance, double omega, int is_aggregated[3])
{
    int nx = dm->nx, ny = dm->ny, nz = dm->nz;
    int i, j, k, iter, ista, offset;
    int is_same[3];
    double rsd0tol = 0.0, rsd_norm = 0.0, temp;
    double *rsd;

    rsd = calloc((nx+2)*(ny+2)*(nz+2), sizeof(double));

    // row-major index
    #define IDX(i,j,k) ((i)*(ny+2)*(nz+2) + (j)*(nz+2) + (k))

    mv_mul_poisson_matrix(rsd, a_poisson, sol, dm, is_aggregated);
    
    #pragma omp parallel for private(i,j,k)
    for(i=1; i<=nx; i++)
        for(j=1; j<=ny; j++)
            for(k=1; k<=nz; k++)
                rsd[IDX(i,j,k)] = rhs[IDX(i,j,k)] - rsd[IDX(i,j,k)];

    vv_dot_3d_matrix(&rsd0tol, rsd, rsd, nx, ny, nz, is_aggregated);
    rsd_norm = rsd0tol;

    // offset = ( (nx%2)*(comm_1d_x.myrank*is_same[0]%2)
    //          + (ny%2)*(comm_1d_y.myrank*is_same[1]%2)
    //          + (nz%2)*(comm_1d_z.myrank*is_same[2]%2) ) % 2;
    offset = ( (nx%2)*(comm_1d_x.myrank%2)
             + (ny%2)*(comm_1d_y.myrank%2)
             + (nz%2)*(comm_1d_z.myrank%2) ) % 2;

    for(iter=0; iter<maxiteration; iter++)
    {
        if ((iter % 10 == 0) && (myrank == 0)) {
            printf("   [RBGS solver] mse: %e, r_mse: %e, rsd0: %e at %d iterations.\n", sqrt(rsd_norm), sqrt(rsd_norm/rsd0tol), sqrt(rsd0tol), iter);
        }
        geometry_halocell_update_selectively(sol, dm, is_aggregated);

        // Red-Black sweep 1
        #pragma omp parallel for private(i,j,k,temp,ista)
        
        for(i=1; i<=nx; i++)
        {
            for(j=1; j<=ny; j++){
                ista = 1 + (i + j + offset)%2;
                for(k=ista; k<=nz; k+=2){
                    temp = *COEFF(a_poisson, 1, i, j, k) * sol[IDX(i-1,j,k)]
                         + *COEFF(a_poisson, 2, i, j, k) * sol[IDX(i+1,j,k)]
                         + *COEFF(a_poisson, 3, i, j, k) * sol[IDX(i,j-1,k)]
                         + *COEFF(a_poisson, 4, i, j, k) * sol[IDX(i,j+1,k)]
                         + *COEFF(a_poisson, 5, i, j, k) * sol[IDX(i,j,k-1)]
                         + *COEFF(a_poisson, 6, i, j, k) * sol[IDX(i,j,k+1)];
                    sol[IDX(i,j,k)] = omega * (rhs[IDX(i,j,k)] - temp) / *COEFF(a_poisson, 0, i, j, k) + (1.0-omega)*sol[IDX(i,j,k)];
                }
            }
        }
        geometry_halocell_update_selectively(sol, dm, is_aggregated);

        // Red-Black sweep 2
        #pragma omp parallel for private(i,j,k,temp,ista)
        
        for(i=1; i<=nx; i++){
            for(j=1; j<=ny; j++){
                ista = 1 + (i + j + 1 + offset)%2;
                for(k=ista; k<=nz; k+=2){
                    temp = *COEFF(a_poisson, 1, i, j, k) * sol[IDX(i-1,j,k)]
                         + *COEFF(a_poisson, 2, i, j, k) * sol[IDX(i+1,j,k)]
                         + *COEFF(a_poisson, 3, i, j, k) * sol[IDX(i,j-1,k)]
                         + *COEFF(a_poisson, 4, i, j, k) * sol[IDX(i,j+1,k)]
                         + *COEFF(a_poisson, 5, i, j, k) * sol[IDX(i,j,k-1)]
                         + *COEFF(a_poisson, 6, i, j, k) * sol[IDX(i,j,k+1)];
                    sol[IDX(i,j,k)] = omega * (rhs[IDX(i,j,k)] - temp) / *COEFF(a_poisson, 0, i, j, k) + (1.0-omega)*sol[IDX(i,j,k)];
                }
            }
        }
   
    // compute final residual
    mv_mul_poisson_matrix(rsd, a_poisson, sol, dm, is_aggregated);

    #pragma omp parallel for private(i,j,k)
    for(i=1; i<=nx; i++)
        for(j=1; j<=ny; j++)
            for(k=1; k<=nz; k++)
                rsd[IDX(i,j,k)] = rhs[IDX(i,j,k)] - rsd[IDX(i,j,k)];

    vv_dot_3d_matrix(&rsd_norm, rsd, rsd, nx, ny, nz, is_aggregated);
    if (sqrt(rsd_norm/rsd0tol) <= tolerance) break;
}
    if(myrank==0)
        printf("[RBGS solver] RBGS solver completed. Iteration = %d, mse = %e, r_mse = %e\n", iter, sqrt(rsd_norm), sqrt(rsd_norm/rsd0tol));
    
    free(rsd);
}

void rbgs_iterator_poisson_matrix(double *sol, matrix_poisson *a_poisson, double *rhs, subdomain *dm, int maxiteration, double omega, int is_aggregated[3]) 
{
    int nx = dm->nx, ny = dm->ny, nz = dm->nz;
    int i, j, k, iter, ista, offset;
    int is_same[3];
    double rsd0tol = 0.0, rsd_norm = 0.0, temp;
    
    double *rsd = (double*)malloc((nx+2)*(ny+2)*(nz+2)*sizeof(double));

    // row-major index
    #define IDX(i,j,k) ((i)*(ny+2)*(nz+2) + (j)*(nz+2) + (k))


    // 1. compute initial residual
    mv_mul_poisson_matrix(rsd, a_poisson, sol, dm, is_aggregated);

    #pragma omp parallel for private(i,j,k)
    for(i=1; i<=nx; i++)
    {
        for(j=1; j<=ny; j++)
        {
            for(k=1; k<=nz; k++)
            {
                rsd[IDX(i,j,k)] = rhs[IDX(i,j,k)] - rsd[IDX(i,j,k)];
            }      
        }
    }
        

    vv_dot_3d_matrix(&rsd0tol, rsd, rsd, nx, ny, nz, is_aggregated);
    rsd_norm = rsd0tol;

    is_same[0] = is_aggregated[0] ? 0 : 1;
    is_same[1] = is_aggregated[1] ? 0 : 1;
    is_same[2] = is_aggregated[2] ? 0 : 1;

    offset = ( (nx%2)*(comm_1d_x.myrank*is_same[0]%2)
             + (ny%2)*(comm_1d_y.myrank*is_same[1]%2)
             + (nz%2)*(comm_1d_z.myrank*is_same[2]%2) ) % 2;


    for(iter=0; iter<maxiteration; iter++){
        geometry_halocell_update_selectively(sol, dm, is_aggregated);
        
        // Red-Black sweep 1
        #pragma omp parallel for private(i,j,k,temp,ista)
        for(i=1; i<=nx; i++){
            for(j=1; j<=ny; j++){
                ista = 1 + (i + j + offset)%2;
                for(k=ista; k<=nz; k+=2){
                    temp = *COEFF(a_poisson, 1, i, j, k) * sol[IDX(i-1,j,k)]
                         + *COEFF(a_poisson, 2, i, j, k) * sol[IDX(i+1,j,k)]
                         + *COEFF(a_poisson, 3, i, j, k) * sol[IDX(i,j-1,k)]
                         + *COEFF(a_poisson, 4, i, j, k) * sol[IDX(i,j+1,k)]
                         + *COEFF(a_poisson, 5, i, j, k) * sol[IDX(i,j,k-1)]
                         + *COEFF(a_poisson, 6, i, j, k) * sol[IDX(i,j,k+1)];
                    sol[IDX(i,j,k)] = omega * (rhs[IDX(i,j,k)] - temp) / *COEFF(a_poisson, 0, i, j, k) + (1.0-omega)*sol[IDX(i,j,k)];
                }
            }
        }

        geometry_halocell_update_selectively(sol, dm, is_aggregated);

        // Red-Black sweep 2
        #pragma omp parallel for private(i,j,k,temp,ista)
        for(i=1; i<=nx; i++){
            for(j=1; j<=ny; j++){
                ista = 1 + (i + j + 1 + offset)%2;
                for(k=ista; k<=nz; k+=2){
                    temp = *COEFF(a_poisson, 1, i, j, k) * sol[IDX(i-1,j,k)]
                         + *COEFF(a_poisson, 2, i, j, k) * sol[IDX(i+1,j,k)]
                         + *COEFF(a_poisson, 3, i, j, k) * sol[IDX(i,j-1,k)]
                         + *COEFF(a_poisson, 4, i, j, k) * sol[IDX(i,j+1,k)]
                         + *COEFF(a_poisson, 5, i, j, k) * sol[IDX(i,j,k-1)]
                         + *COEFF(a_poisson, 6, i, j, k) * sol[IDX(i,j,k+1)];
                    sol[IDX(i,j,k)] = omega * (rhs[IDX(i,j,k)] - temp) / *COEFF(a_poisson, 0, i, j, k) + (1.0-omega)*sol[IDX(i,j,k)];
                }
            }
        }
    }

    // compute final residual
    mv_mul_poisson_matrix(rsd, a_poisson, sol, dm, is_aggregated);

    #pragma omp parallel for private(i,j,k)
    for(i=1; i<=nx; i++)
        for(j=1; j<=ny; j++)
            for(k=1; k<=nz; k++)
                rsd[IDX(i,j,k)] = rhs[IDX(i,j,k)] - rsd[IDX(i,j,k)];

    vv_dot_3d_matrix(&rsd_norm, rsd, rsd, nx, ny, nz, is_aggregated);

    if(myrank==0)
        printf("[RBGS iterator] RBGS iterator completed. Iteration = %d, mse = %e, r_mse = %e, rsd0 = %e \n", iter, sqrt(rsd_norm), sqrt(rsd_norm/rsd0tol), sqrt(rsd0tol));

    free(rsd);
}