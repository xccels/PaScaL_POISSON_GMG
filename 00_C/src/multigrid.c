#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "multigrid.h"
// #include "mpi_topology.h"   
// #include "matrix.h" 
// #include "multigrid_common.h" 
#include "poisson_matrix_operator.h" 
// #include "rbgs_poisson_matrix.h"     

static int lv_gdm_coarsest_max;
static int lv_gdm_coarsest_x;
static int lv_gdm_coarsest_y;
static int lv_gdm_coarsest_z;
static int n_levels;
static int n_vcycles;
static int lv_aggregation;
static int lv_aggregation_x;
static int lv_aggregation_y;
static int lv_aggregation_z;
static int lv_aggregation_max;
static int aggregation_type;

static subdomain *mg_sdm = NULL;         
static matrix_poisson *mg_a_poisson = NULL;

#define IDX(i,j,k, nx,ny,nz) ((i)*(ny+2)*(nz+2) + (j)*(nz+2) + (k))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

void multigrid_create(subdomain *sdm, int nlevel, int ncycle, int aggr_method, int aggr_level)
{

    int l, ierr;
    
    n_levels = nlevel;
    n_vcycles = ncycle;
    lv_aggregation = aggr_level;
    aggregation_type = aggr_method;

    if (n_levels <= 1) {
        if (myrank == 0) {
            printf("[Error] The number of levels should be larger than 1. Current number: %d\n", n_levels);
        }
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    mg_sdm = (subdomain *)malloc((n_levels+1) * sizeof(subdomain));
    if (!mg_sdm) {
        if (myrank == 0) printf("[Error] Failed to allocate mg_sdm\n");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    switch (aggregation_type) {
        case 0:
            multigrid_subdomain_create_no_aggregation(sdm);
            break;
        case 1:
            multigrid_subdomain_create_single_aggregation(sdm);
            break;
        case 2:
            multigrid_subdomain_create_adaptive_aggregation(sdm);
            break;
        default:
            if (myrank == 0) printf("[Error] Aggregation method should be 0, 1, or 2. Current: %d\n", aggregation_type);
            MPI_Finalize();
            exit(EXIT_FAILURE);
    }

    multigrid_allocate_subdomain_variables();

    if (myrank == 0) {
        printf("[MG] Aggregation info. of subdomain: %c %c %c\n",
           sdm->is_aggregated[0] ? 'T' : 'F',
           sdm->is_aggregated[1] ? 'T' : 'F',
           sdm->is_aggregated[2] ? 'T' : 'F');
    }

    for (l = 1; l <= n_levels; l++) {
        if (myrank == 0) {
            printf("[MG] Aggregation info. of level %d: %c %c %c\n",l, 
                mg_sdm[l].is_aggregated[0] ? 'T' : 'F',
                mg_sdm[l].is_aggregated[1] ? 'T' : 'F',
                mg_sdm[l].is_aggregated[2] ? 'T' : 'F');
        }
    }

    if (myrank == 0) printf("[MG] Multigrid geometry constructed.\n");

    mg_a_poisson = (matrix_poisson *)malloc((n_levels+1) * sizeof(matrix_poisson));
    if (!mg_a_poisson) {
        if (myrank == 0) printf("[Error] Failed to allocate mg_a_poisson\n");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    for (l = 1; l <= n_levels; l++) {
        matrix_poisson_create(&mg_a_poisson[l], &mg_sdm[l]);
    }

    if (myrank == 0) {printf("[MG] Poisson matrix in multigrid constructed.\n");}


#ifdef DEBUG_GRID
    multigrid_common_print_subdomain_grid_info(sdm, mg_sdm, n_levels);
#endif

#ifdef DEBUG_MATRIX
    for (l = 1; l <= n_levels; l++) {
        multigrid_common_print_poisson_matrix_info(&mg_sdm[l], &mg_a_poisson[l], l);
    }
#endif

}

void multigrid_subdomain_create_no_aggregation(subdomain *sdm) 
{
    int l;
    int nx, ny, nz;
    double *dxm_f, *dxg_f, *xg_f;
    double *dym_f, *dyg_f, *yg_f;
    double *dzm_f, *dzg_f, *zg_f;

    if (myrank == 0) {printf("[MG] Grid coarsening without aggregation.\n");}

    if (lv_aggregation != 0) {
        if (myrank == 0) printf("[Error] Aggregation level should be 0 for no aggregation method.\n");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    nx = sdm->nx;
    ny = sdm->ny;
    nz = sdm->nz;

    // aggregation check for processes
    for (l = 1; l <= n_levels; l++) {
        mg_sdm[l].is_aggregated[0] = (comm_1d_x.nprocs == 1);
        mg_sdm[l].is_aggregated[1] = (comm_1d_y.nprocs == 1);
        mg_sdm[l].is_aggregated[2] = (comm_1d_z.nprocs == 1);
    }

    for (l = 1; l <= n_levels; l++) {
        if (myrank == 0) printf("[MG] Grid coarsening at level %d\n", l);

        // X-direction
        if (nx % 2 == 0) {
            nx = nx/2;
            lv_gdm_coarsest_x = l;
            if (myrank == 0)
                printf("[MG] X-grid coarsening at level %d: %d grids reduced to %d grids.\n",
                       l, nx, nx/2);    // nx到nx/2 还是 nx*2到nx？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？
        } else {
            if (mg_sdm[l].is_aggregated[0]) {
                if (myrank == 0) {
                    printf("[MG] No more X-grid coarsening from level %d. Keeping X-grid number.\n", l);
                    printf("[MG] The number of grid is %d and number processes is %d\n",
                           nx, comm_1d_x.nprocs);
                }
            } else {
                if (myrank == 0) {
                    printf("[Error] X-grid coarsening impossible for multiple processes from level %d\n", l);
                    printf("[Error] The number of grid per process is %d and number processes is %d\n",
                           nx, comm_1d_x.nprocs);
                }
                MPI_Finalize();
                exit(EXIT_FAILURE);
            }
        }

        // Y-direction
        if (ny % 2 == 0) {
            ny = ny/2;
            lv_gdm_coarsest_y = l;
            if (myrank == 0)
                printf("[MG] Y-grid coarsening at level %d: %d grids reduced to %d grids.\n",
                       l, ny, ny/2);    // ny到ny/2 还是 ny*2到ny？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？
        } else {
            if (mg_sdm[l].is_aggregated[1]) {
                if (myrank == 0) {
                    printf("[MG] No more Y-grid coarsening from level %d. Keeping Y-grid number.\n", l);
                    printf("[MG] The number of grid is %d and number processes is %d\n",
                           ny, comm_1d_y.nprocs);
                }
            } else {
                if (myrank == 0) {
                    printf("[Error] Y-grid coarsening impossible for multiple processes from level %d\n", l);
                    printf("[Error] The number of grid per process is %d and number processes is %d\n",
                           ny, comm_1d_y.nprocs);
                }
                MPI_Finalize();
                exit(EXIT_FAILURE);
            }
        }

        // Z-direction
        if (nz % 2 == 0) {
            nz = nz/2;
            lv_gdm_coarsest_z = l;
            if (myrank == 0)
                printf("[MG] Z-grid coarsening at level %d: %d grids reduced to %d grids.\n",
                       l, nz, nz/2);    // nz到nz/2 还是 nz*2到nz？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？
        } else {
            if (mg_sdm[l].is_aggregated[2]) {
                if (myrank == 0) {
                    printf("[MG] No more Z-grid coarsening from level %d. Keeping Z-grid number.\n", l);
                    printf("[MG] The number of grid is %d and number processes is %d\n",
                           nz, comm_1d_z.nprocs);
                }
            } else {
                if (myrank == 0) {
                    printf("[Error] Z-grid coarsening impossible for multiple processes from level %d\n", l);
                    printf("[Error] The number of grid per process is %d and number processes is %d\n",
                           nz, comm_1d_z.nprocs);
                }
                MPI_Finalize();
                exit(EXIT_FAILURE);
            }
        }

        mg_sdm[l].nx = nx;
        mg_sdm[l].ny = ny;
        mg_sdm[l].nz = nz;
    }

    lv_gdm_coarsest_max = MAX(MAX(lv_gdm_coarsest_x, lv_gdm_coarsest_y), lv_gdm_coarsest_z);
    if (myrank == 0) {
        printf("[MG] Number of levels : %d\n", n_levels);
        printf("[MG] Final coarsest levels max, x, y, z directions : %d %d %d %d\n",
               lv_gdm_coarsest_max, lv_gdm_coarsest_x, lv_gdm_coarsest_y, lv_gdm_coarsest_z);
    }

    // Check the number of levels and the max coarsest level
    if (n_levels > lv_gdm_coarsest_max) {
        if (myrank == 0) {
            printf("[MG] Number of levels is greater than the coarsest level in global domain.\n");
            printf("[MG] It is not allowed and reduce the number of levels.\n");
        }
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    if (myrank == 0) printf("[MG] Generating grid dimension\n");

    // In x-direction
    dxm_f = sdm->dxm;
    dxg_f = sdm->dxg;
    xg_f  = sdm->xg;

    for (l = 1; l <= n_levels; l++) {
        nx = mg_sdm[l].nx;
        mg_sdm[l].dxm = (double*)malloc((nx+2) * sizeof(double));
        mg_sdm[l].dxg = (double*)malloc((nx+2) * sizeof(double));
        mg_sdm[l].xg  = (double*)malloc((nx+2) * sizeof(double));

        multigrid_subdomain_no_aggregation_make_grid(mg_sdm[l].dxm, mg_sdm[l].dxg, mg_sdm[l].xg, nx, dxm_f, dxg_f, xg_f, sdm->ox, l, lv_gdm_coarsest_x, comm_1d_x, 'x');

        dxm_f = NULL;
        dxg_f = NULL;
        xg_f  = NULL;

        dxm_f = mg_sdm[l].dxm;
        dxg_f = mg_sdm[l].dxg;
        xg_f  = mg_sdm[l].xg;
    }
    dxm_f = NULL;
    dxg_f = NULL;
    xg_f  = NULL;

    // In y-direction
    dym_f = sdm->dym;
    dyg_f = sdm->dyg;
    yg_f  = sdm->yg;

    for (l = 1; l <= n_levels; l++) {
        ny = mg_sdm[l].ny;
        mg_sdm[l].dym = (double*)malloc((ny+2) * sizeof(double));
        mg_sdm[l].dyg = (double*)malloc((ny+2) * sizeof(double));
        mg_sdm[l].yg  = (double*)malloc((ny+2) * sizeof(double));

        multigrid_subdomain_no_aggregation_make_grid(mg_sdm[l].dym, mg_sdm[l].dyg, mg_sdm[l].yg, ny, dym_f, dyg_f, yg_f, sdm->oy, l, lv_gdm_coarsest_y, comm_1d_y, 'y');

        dym_f = NULL;
        dyg_f = NULL;
        yg_f  = NULL;
        
        dym_f = mg_sdm[l].dym;
        dyg_f = mg_sdm[l].dyg;
        yg_f  = mg_sdm[l].yg;
    }
    dym_f = NULL;
    dyg_f = NULL;
    yg_f  = NULL;

    // In z-direction
    dzm_f = sdm->dzm;
    dzg_f = sdm->dzg;
    zg_f  = sdm->zg;

    for (l = 1; l <= n_levels; l++) {
        nz = mg_sdm[l].nz;
        mg_sdm[l].dzm = (double*)malloc((nz+2) * sizeof(double));
        mg_sdm[l].dzg = (double*)malloc((nz+2) * sizeof(double));
        mg_sdm[l].zg  = (double*)malloc((nz+2) * sizeof(double));

        multigrid_subdomain_no_aggregation_make_grid(mg_sdm[l].dzm, mg_sdm[l].dzg, mg_sdm[l].zg, nz, dzm_f, dzg_f, zg_f, sdm->oz, l, lv_gdm_coarsest_z, comm_1d_z, 'z');

        dzm_f = NULL;
        dzg_f = NULL;
        zg_f  = NULL;
        
        dzm_f = mg_sdm[l].dzm;
        dzg_f = mg_sdm[l].dzg;
        zg_f  = mg_sdm[l].zg;
    }
    dzm_f = NULL;
    dzg_f = NULL;
    zg_f  = NULL;

}

void multigrid_subdomain_create_single_aggregation(subdomain *sdm) 
{
    int l;
    int nx, ny, nz;
    int nx_lv_aggregation, ny_lv_aggregation, nz_lv_aggregation;

    double *dxm_f, *dxg_f, *xg_f;
    double *dym_f, *dyg_f, *yg_f;
    double *dzm_f, *dzg_f, *zg_f;

    if (myrank == 0) {printf("[MG] Grid coarsening with aggretation level = %d\n", lv_aggregation);}

    lv_aggregation_x = lv_aggregation;
    lv_aggregation_y = lv_aggregation;
    lv_aggregation_z = lv_aggregation;
    
    if (lv_aggregation == 0) {
        if (myrank == 0) printf("[Error] Aggretation level should be larger than 0.\n");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    nx = sdm->nx;
    ny = sdm->ny;
    nz = sdm->nz;


    for (l = 1; l <= n_levels; l++) {
        if (myrank == 0) printf("[MG] Grid coarsening at level %d\n", l);

        if (l < lv_aggregation) 
        {
            // X
            if (nx % 2 == 1) 
            {
                if (myrank == 0) 
                {
                    printf("[Error] X-grid coarsening impossible at level = %d, less than aggregation level = %d\n", l, lv_aggregation);
                    printf("[Error] The number of grid per process is %d and number processes is %d\n", nx, comm_1d_x.nprocs);
                }
                MPI_Finalize();
                exit(1);
            } 
            else {
                mg_sdm[l].nx = nx / 2;
                if (myrank == 0)
                {
                    printf("[MG] X-grid coarsening at level %d: %d grids reduced to %d grids.\n", l, nx, nx / 2);
                }
            }
            // Y
            if (ny % 2 == 1) 
            {
                if (myrank == 0) {
                    printf("[Error] Y-grid coarsening impossible at level = %d, less than aggregation level = %d\n", l, lv_aggregation);
                    printf("[Error] The number of grid per process is %d and number processes is %d\n", ny, comm_1d_y.nprocs);
                }
                MPI_Finalize();
                exit(1);
            } 
            else {
                mg_sdm[l].ny = ny / 2;
                if (myrank == 0)
                {
                    printf("[MG] Y-grid coarsening at level %d: %d grids reduced to %d grids.\n", l, ny, ny / 2);
                }
            }
            // Z
            if (nz % 2 == 1) {
                if (myrank == 0) {
                    printf("[Error] Z-grid coarsening impossible at level = %d, less than aggregation level = %d\n", l, lv_aggregation);
                    printf("[Error] The number of grid per process is %d and number processes is %d\n", nz, comm_1d_z.nprocs);
                }
                MPI_Finalize();
                exit(1);
            } 
            else {
                mg_sdm[l].nz = nz / 2;
                if (myrank == 0)
                {
                    printf("[MG] Z-grid coarsening at level %d: %d grids reduced to %d grids.\n", l, nz, nz / 2);
                }
            }
        }
        else if (l == lv_aggregation)
        {
            // X
            if (nx % 2 == 1) 
            {
                if (myrank == 0) 
                {
                    printf("[Error] X-grid coarsening impossible at level = %d, equal to aggregation level = %d\n", l, lv_aggregation);
                    printf("[Error] The number of grid per process is %d and number processes is %d\n", nx, comm_1d_x.nprocs);
                }
                MPI_Finalize();
                exit(1);
            } 
            else {
                nx_lv_aggregation = nx / 2;
                mg_sdm[l].nx = (nx / 2) * comm_1d_x.nprocs;
                lv_gdm_coarsest_x = l;
                if (myrank == 0)
                {
                    printf("[MG] X-grid coarsening at level %d: %d grids reduced to %d grids.\n", l, nx, nx / 2);
                    printf("[MG] X-grid aggregation at level %d: %d grids aggregated to %d grids.\n", l, nx / 2, mg_sdm[l].nx);
                }
            }
            // Y
            if (ny % 2 == 1) 
            {
                if (myrank == 0) {
                    printf("[Error] Y-grid coarsening impossible at level = %d, equal to aggregation level = %d\n", l, lv_aggregation);
                    printf("[Error] The number of grid per process is %d and number processes is %d\n", ny, comm_1d_y.nprocs);
                }
                MPI_Finalize();
                exit(1);
            } 
            else {
                ny_lv_aggregation = ny / 2;
                mg_sdm[l].ny = (ny / 2) * comm_1d_y.nprocs;
                lv_gdm_coarsest_y = l;
                if (myrank == 0) {
                    printf("[MG] Y-grid coarsening at level %d: %d grids reduced to %d grids.\n", l, ny, ny / 2);
                    printf("[MG] Y-grid aggregation at level %d: %d grids aggregated to %d grids.\n", l, ny / 2, mg_sdm[l].ny);
                }
            }
            // Z
            if (nz % 2 == 1) {
                if (myrank == 0) {
                    printf("[Error] Z-grid coarsening impossible at level = %d, equal to aggregation level = %d\n", l, lv_aggregation);
                    printf("[Error] The number of grid per process is %d and number processes is %d\n", nz, comm_1d_z.nprocs);
                }
                MPI_Finalize();
                exit(1);
            } 
            else {
                nz_lv_aggregation = nz / 2;
                mg_sdm[l].nz = (nz / 2) * comm_1d_z.nprocs;
                lv_gdm_coarsest_z = l;
                if (myrank == 0) {
                    printf("[MG] Z-grid coarsening at level %d: %d grids reduced to %d grids.\n", l, nz, nz / 2);
                    printf("[MG] Z-grid aggregation at level %d: %d grids aggregated to %d grids.\n", l, nz / 2, mg_sdm[l].nz);
                }
            }
        }
        else if (l > lv_aggregation)
        {
            // X
            if (nx % 2 == 1) 
            {
                if (myrank == 0) 
                {
                    printf("[MG] No more X-grid coarsening from level %d. Keeping X-grid number.\n", l);
                    printf("[MG] The number of grid is %d and number processes is %d\n", nx, comm_1d_x.nprocs);
                }
                mg_sdm[l].nx = nx;
            } 
            else {
                if (myrank == 0){
                    printf("[MG] X-grid coarsening at level %d: %d grids reduced to %d grids.\n", l, nx, nx / 2);
                }
                mg_sdm[l].nx = nx / 2;
                lv_gdm_coarsest_x = l;
            }
            // Y
            if (ny % 2 == 1) 
            {
                if (myrank == 0) 
                {
                    printf("[MG] No more Y-grid coarsening from level %d. Keeping Y-grid number.\n", l);
                    printf("[MG] The number of grid is %d and number processes is %d\n", ny, comm_1d_y.nprocs);
                }
                mg_sdm[l].ny = ny;
            } 
            else {
                if (myrank == 0){
                    printf("[MG] Y-grid coarsening at level %d: %d grids reduced to %d grids.\n", l, ny, ny / 2);
                }
                mg_sdm[l].ny = ny / 2;
                lv_gdm_coarsest_y = l;
            }
            // Z
            if (nz % 2 == 1) {
                if (myrank == 0) 
                {
                    printf("[MG] No more Z-grid coarsening from level %d. Keeping Z-grid number.\n", l);
                    printf("[MG] The number of grid is %d and number processes is %d\n", nz, comm_1d_z.nprocs);
                }
                 mg_sdm[l].nz = nz;
            } 
            else {
                if (myrank == 0){
                    printf("[MG] Z-grid coarsening at level %d: %d grids reduced to %d grids.\n", l, nz, nz / 2);
                }
                mg_sdm[l].nz = nz / 2;
                lv_gdm_coarsest_z = l;
            }
        }

        nx = mg_sdm[l].nx;
        ny = mg_sdm[l].ny;
        nz = mg_sdm[l].nz;
        if ((nx*ny*nz)%2 == 1)
        {
            break;
        }
    }

    for (int l = 1; l <= n_levels; l++) 
    {
        if (l < lv_aggregation) 
        {
            mg_sdm[l].is_aggregated[0] = (comm_1d_x.nprocs == 1);
            mg_sdm[l].is_aggregated[1] = (comm_1d_y.nprocs == 1);
            mg_sdm[l].is_aggregated[2] = (comm_1d_z.nprocs == 1);
        } 
        else if (l >= lv_aggregation) 
        {
            mg_sdm[l].is_aggregated[0] = 1;
            mg_sdm[l].is_aggregated[1] = 1;
            mg_sdm[l].is_aggregated[2] = 1;
        }
    }

    lv_gdm_coarsest_max = MAX(MAX(lv_gdm_coarsest_x, lv_gdm_coarsest_y), lv_gdm_coarsest_z);
    if (myrank == 0) {
        printf("[MG] Number of levels : %d\n", n_levels);
        printf("[MG] Final coarsest levels max, x, y, z directions : %d %d %d %d\n",
               lv_gdm_coarsest_max, lv_gdm_coarsest_x, lv_gdm_coarsest_y, lv_gdm_coarsest_z);
    }

    // Check the number of levels and the max coarsest level
    if (n_levels > lv_gdm_coarsest_max) {
        if (myrank == 0) {
            printf("[MG] Number of levels is greater than the coarsest level in global domain.\n");
            printf("[MG] It is not allowed and reduce the number of levels.\n");
        }
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    if (myrank == 0){
        printf("[MG] Final gather level = %d\n",lv_aggregation);
    }
    if (lv_aggregation > n_levels)
    {
        if (myrank == 0){
            printf("[MG] Gather level is greater than the coarsest level and no gather level exists\n");
        }
    }
    if (myrank == 0) printf("[MG] Generating grid dimension\n");

    // In x-direction
    dxm_f = sdm->dxm;
    dxg_f = sdm->dxg;
    xg_f  = sdm->xg;

    for (l = 1; l <= n_levels; l++) {
        nx = mg_sdm[l].nx;
        mg_sdm[l].dxm = (double*)malloc((nx+2) * sizeof(double));
        mg_sdm[l].dxg = (double*)malloc((nx+2) * sizeof(double));
        mg_sdm[l].xg  = (double*)malloc((nx+2) * sizeof(double));

        multigrid_subdomain_aggregation_make_grid(mg_sdm[l].dxm, mg_sdm[l].dxg, mg_sdm[l].xg, nx, dxm_f, dxg_f, xg_f, sdm->ox, l, lv_gdm_coarsest_x, lv_aggregation_x, nx_lv_aggregation, comm_1d_x, 'x');

        dxm_f = NULL;
        dxg_f = NULL;
        xg_f  = NULL;

        dxm_f = mg_sdm[l].dxm;
        dxg_f = mg_sdm[l].dxg;
        xg_f  = mg_sdm[l].xg;
    }
    dxm_f = NULL;
    dxg_f = NULL;
    xg_f  = NULL;

    // In y-direction
    dym_f = sdm->dym;
    dyg_f = sdm->dyg;
    yg_f  = sdm->yg;

    for (l = 1; l <= n_levels; l++) {
        ny = mg_sdm[l].ny;
        mg_sdm[l].dym = (double*)malloc((ny+2) * sizeof(double));
        mg_sdm[l].dyg = (double*)malloc((ny+2) * sizeof(double));
        mg_sdm[l].yg  = (double*)malloc((ny+2) * sizeof(double));

        multigrid_subdomain_aggregation_make_grid(mg_sdm[l].dym, mg_sdm[l].dyg, mg_sdm[l].yg, ny, dym_f, dyg_f, yg_f, sdm->oy, l, lv_gdm_coarsest_y, lv_aggregation_y, ny_lv_aggregation, comm_1d_y, 'y');

        dym_f = NULL;
        dyg_f = NULL;
        yg_f  = NULL;
        
        dym_f = mg_sdm[l].dym;
        dyg_f = mg_sdm[l].dyg;
        yg_f  = mg_sdm[l].yg;
    }
    dym_f = NULL;
    dyg_f = NULL;
    yg_f  = NULL;

    // In z-direction
    dzm_f = sdm->dzm;
    dzg_f = sdm->dzg;
    zg_f  = sdm->zg;

    for (l = 1; l <= n_levels; l++) {
        nz = mg_sdm[l].nz;
        mg_sdm[l].dzm = (double*)malloc((nz+2) * sizeof(double));
        mg_sdm[l].dzg = (double*)malloc((nz+2) * sizeof(double));
        mg_sdm[l].zg  = (double*)malloc((nz+2) * sizeof(double));

        multigrid_subdomain_aggregation_make_grid(mg_sdm[l].dzm, mg_sdm[l].dzg, mg_sdm[l].zg, nz, dzm_f, dzg_f, zg_f, sdm->oz, l, lv_gdm_coarsest_z, lv_aggregation_z, nz_lv_aggregation, comm_1d_z, 'z');

        dzm_f = NULL;
        dzg_f = NULL;
        zg_f  = NULL;
        
        dzm_f = mg_sdm[l].dzm;
        dzg_f = mg_sdm[l].dzg;
        zg_f  = mg_sdm[l].zg;
    }
    dzm_f = NULL;
    dzg_f = NULL;
    zg_f  = NULL;


}

void multigrid_subdomain_create_adaptive_aggregation(subdomain *sdm) 
{
    int l;
    int nx, ny, nz;
    int nx_lv_aggregation, ny_lv_aggregation, nz_lv_aggregation;

    double *dxm_f, *dxg_f, *xg_f;
    double *dym_f, *dyg_f, *yg_f;
    double *dzm_f, *dzg_f, *zg_f;

    if (myrank == 0) {printf("[MG] Grid coarsening with adaptive aggretation. lv_aggregation is neglected\n");}

    nx = sdm->nx;
    ny = sdm->ny;
    nz = sdm->nz;

    lv_aggregation_x = n_levels+1;
    lv_aggregation_y = n_levels+1;
    lv_aggregation_z = n_levels+1;


    for (l = 1; l <= n_levels; l++) {
        if (myrank == 0) printf("[MG] Grid coarsening at level %d\n", l);
    
        // X
        if (nx % 2 == 1) 
        {
            mg_sdm[l].nx = nx;
            if (myrank == 0) 
            {
                printf("[MG] No more x-grid coarsening at level %d in serial. Keeping x-grid number\n", l);
                printf("[MG] Aggregation level in x-direction is %d\n", lv_aggregation_x);
                printf("[MG] The number of grid is %d and number processes is %d\n", mg_sdm[l].nx, comm_1d_x.nprocs);
            }
        } 
        else {
            if( ((nx/2) % 2 == 1) && (lv_aggregation_x == (n_levels+1)) )
            {
                if (comm_1d_x.nprocs == 1)
                {
                    lv_aggregation_x = 0;
                }
                else{
                    lv_aggregation_x = l;
                }
                lv_gdm_coarsest_x = l;
                nx_lv_aggregation = nx/2;
                mg_sdm[l].nx = nx/2 * comm_1d_x.nprocs;
                if (myrank == 0)
                {
                    printf("[MG] x-grid is odd number at level %d in parallel.\n", l);
                    printf("[MG] Aggregation level in x-direction prescribed as %d\n", lv_aggregation_x);
                    printf("[MG] The number of grid is %d and number processes is %d\n", mg_sdm[l].nx, comm_1d_x.nprocs);
                }
            }
            else{
                mg_sdm[l].nx = nx/2;
                lv_gdm_coarsest_x = l;
                if (myrank == 0) 
                {
                    printf("[MG] x-grid coarsening at level %d: %d grids reduced to %d grids.\n", l, nx, mg_sdm[l].nx);
                }
            }  
        }

        // Y
        if (ny % 2 == 1) 
        {
            mg_sdm[l].ny = ny;
            if (myrank == 0) 
            {
                printf("[MG] No more y-grid coarsening at level %d in serial. Keeping y-grid number\n", l);
                printf("[MG] Aggregation level in y-direction is %d\n", lv_aggregation_y);
                printf("[MG] The number of grid is %d and number processes is %d\n", mg_sdm[l].ny, comm_1d_y.nprocs);
            }
        } 
        else {
            if( ((ny/2) % 2 == 1) && (lv_aggregation_y == (n_levels+1)) )
            {
                if (comm_1d_y.nprocs == 1)
                {
                    lv_aggregation_y = 0;
                }
                else{
                    lv_aggregation_y = l;
                }
                lv_gdm_coarsest_y = l;
                ny_lv_aggregation = ny/2;
                mg_sdm[l].ny = ny/2 * comm_1d_y.nprocs;
                if (myrank == 0)
                {
                    printf("[MG] y-grid is odd number at level %d in parallel.\n", l);
                    printf("[MG] Aggregation level in y-direction prescribed as %d\n", lv_aggregation_y);
                    printf("[MG] The number of grid is %d and number processes is %d\n", mg_sdm[l].ny, comm_1d_y.nprocs);
                }
            }
            else{
                mg_sdm[l].ny = ny/2;
                lv_gdm_coarsest_y = l;
                if (myrank == 0) 
                {
                    printf("[MG] y-grid coarsening at level %d: %d grids reduced to %d grids.\n", l, ny, mg_sdm[l].ny);
                }
            }  
        }

        // Z
        if (nz % 2 == 1) 
        {
            mg_sdm[l].nz = nz;
            if (myrank == 0) 
            {
                printf("[MG] No more z-grid coarsening at level %d in serial. Keeping z-grid number\n", l);
                printf("[MG] Aggregation level in z-direction is %d\n", lv_aggregation_z);
                printf("[MG] The number of grid is %d and number processes is %d\n", mg_sdm[l].nz, comm_1d_z.nprocs);
            }
        } 
        else {
            if( ((nz/2) % 2 == 1) && (lv_aggregation_z == (n_levels+1)) )
            {
                if (comm_1d_z.nprocs == 1)
                {
                    lv_aggregation_z = 0;
                }
                else{
                    lv_aggregation_z = l;
                }
                lv_gdm_coarsest_z = l;
                nz_lv_aggregation = nz/2;
                mg_sdm[l].nz = nz/2 * comm_1d_z.nprocs;
                if (myrank == 0)
                {
                    printf("[MG] z-grid is odd number at level %d in parallel.\n", l);
                    printf("[MG] Aggregation level in z-direction prescribed as %d\n", lv_aggregation_z);
                    printf("[MG] The number of grid is %d and number processes is %d\n", mg_sdm[l].nz, comm_1d_z.nprocs);
                }
            }
            else{
                mg_sdm[l].nz = nz/2;
                lv_gdm_coarsest_z = l;
                if (myrank == 0) 
                {
                    printf("[MG] z-grid coarsening at level %d: %d grids reduced to %d grids.\n", l, nz, mg_sdm[l].nz);
                }
            }  
        }

        nx = mg_sdm[l].nx;
        ny = mg_sdm[l].ny;
        nz = mg_sdm[l].nz;
        if ((nx*ny*nz)%2 == 1)
        {
            break;
        }
    }


    for (int l = 1; l <= n_levels; l++) 
    {
        if (l < lv_aggregation_x) 
        {
            mg_sdm[l].is_aggregated[0] = (comm_1d_x.nprocs == 1);
        } 
        else if (l >= lv_aggregation_x) 
        {
            mg_sdm[l].is_aggregated[0] = 1;
        }

        if (l < lv_aggregation_y) 
        {
            mg_sdm[l].is_aggregated[1] = (comm_1d_y.nprocs == 1);
        } 
        else if (l >= lv_aggregation_y) 
        {
            mg_sdm[l].is_aggregated[1] = 1;
        }

        if (l < lv_aggregation_z) 
        {
            mg_sdm[l].is_aggregated[2] = (comm_1d_z.nprocs == 1);
        } 
        else if (l >= lv_aggregation_z) 
        {
            mg_sdm[l].is_aggregated[2] = 1;
        }
    }

    lv_gdm_coarsest_max = MAX(MAX(lv_gdm_coarsest_x, lv_gdm_coarsest_y), lv_gdm_coarsest_z);
    lv_aggregation_max = MAX(MAX(lv_aggregation_x, lv_aggregation_y), lv_aggregation_z);
    if (myrank == 0) {
        printf("[MG] Number of levels : %d\n", n_levels);
        printf("[MG] Final aggregation levels max, x, y, z directions : %d %d %d %d\n",
               lv_aggregation_max, lv_aggregation_x, lv_aggregation_y, lv_aggregation_z);
        printf("[MG] Final coarsest levels max, x, y, z directions : %d %d %d %d\n",
               lv_gdm_coarsest_max, lv_gdm_coarsest_x, lv_gdm_coarsest_y, lv_gdm_coarsest_z);
    }

    // Check the number of levels and the max coarsest level
    if (n_levels != lv_gdm_coarsest_max) {
        if (myrank == 0) {
            printf("[MG] Number of levels should be equal to the coarsest level in global domain.\n");
            printf("[MG] It is not allowed and reduce the number of levels.\n");
        }
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }


    // In x-direction
    dxm_f = sdm->dxm;
    dxg_f = sdm->dxg;
    xg_f  = sdm->xg;

    for (l = 1; l <= n_levels; l++) {
        nx = mg_sdm[l].nx;
        mg_sdm[l].dxm = (double*)malloc((nx+2) * sizeof(double));
        mg_sdm[l].dxg = (double*)malloc((nx+2) * sizeof(double));
        mg_sdm[l].xg  = (double*)malloc((nx+2) * sizeof(double));

        multigrid_subdomain_aggregation_make_grid(mg_sdm[l].dxm, mg_sdm[l].dxg, mg_sdm[l].xg, nx, dxm_f, dxg_f, xg_f, sdm->ox, l, lv_gdm_coarsest_x, lv_aggregation_x, nx_lv_aggregation, comm_1d_x, 'x');

        dxm_f = NULL;
        dxg_f = NULL;
        xg_f  = NULL;

        dxm_f = mg_sdm[l].dxm;
        dxg_f = mg_sdm[l].dxg;
        xg_f  = mg_sdm[l].xg;
    }
    dxm_f = NULL;
    dxg_f = NULL;
    xg_f  = NULL;

    // In y-direction
    dym_f = sdm->dym;
    dyg_f = sdm->dyg;
    yg_f  = sdm->yg;

    for (l = 1; l <= n_levels; l++) {
        ny = mg_sdm[l].ny;
        mg_sdm[l].dym = (double*)malloc((ny+2) * sizeof(double));
        mg_sdm[l].dyg = (double*)malloc((ny+2) * sizeof(double));
        mg_sdm[l].yg  = (double*)malloc((ny+2) * sizeof(double));

        multigrid_subdomain_aggregation_make_grid(mg_sdm[l].dym, mg_sdm[l].dyg, mg_sdm[l].yg, ny, dym_f, dyg_f, yg_f, sdm->oy, l, lv_gdm_coarsest_y, lv_aggregation_y, ny_lv_aggregation, comm_1d_y, 'y');

        dym_f = NULL;
        dyg_f = NULL;
        yg_f  = NULL;
        
        dym_f = mg_sdm[l].dym;
        dyg_f = mg_sdm[l].dyg;
        yg_f  = mg_sdm[l].yg;
    }
    dym_f = NULL;
    dyg_f = NULL;
    yg_f  = NULL;

    // In z-direction
    dzm_f = sdm->dzm;
    dzg_f = sdm->dzg;
    zg_f  = sdm->zg;

    for (l = 1; l <= n_levels; l++) {
        nz = mg_sdm[l].nz;
        mg_sdm[l].dzm = (double*)malloc((nz+2) * sizeof(double));
        mg_sdm[l].dzg = (double*)malloc((nz+2) * sizeof(double));
        mg_sdm[l].zg  = (double*)malloc((nz+2) * sizeof(double));

        multigrid_subdomain_aggregation_make_grid(mg_sdm[l].dzm, mg_sdm[l].dzg, mg_sdm[l].zg, nz, dzm_f, dzg_f, zg_f, sdm->oz, l, lv_gdm_coarsest_z, lv_aggregation_z, nz_lv_aggregation, comm_1d_z, 'z');

        dzm_f = NULL;
        dzg_f = NULL;
        zg_f  = NULL;
        
        dzm_f = mg_sdm[l].dzm;
        dzg_f = mg_sdm[l].dzg;
        zg_f  = mg_sdm[l].zg;
    }
    dzm_f = NULL;
    dzg_f = NULL;
    zg_f  = NULL;

}


void multigrid_subdomain_no_aggregation_make_grid(double *dxm, double *dxg, double *xg, int nx, const double *dxm_f, const double *dxg_f, const double *xg_f, double ox, int lv_cur, int lv_coarsest, cart_comm_1d comm_1d, char dir)
{
    int i;
    MPI_Request request1[2], request2[2];

    if (myrank == 0) {
        printf("[MG] level : %d, dir : %c\n", lv_cur, dir);
    }

    for (i = 0; i <= nx+1; i++) {
        dxm[i] = 0.0;
        dxg[i] = 0.0;
        xg[i]  = 0.0;
    }

    if (lv_cur <= lv_coarsest) {
        for (i = 1; i <= nx; i++) {
            dxm[i] = dxm_f[2*i-1] + dxm_f[2*i];
        }

        if (comm_1d.west_rank != MPI_PROC_NULL) {
            MPI_Isend(&dxm[1], 1, MPI_DOUBLE, comm_1d.west_rank, 1, comm_1d.mpi_comm, &request1[0]);
            MPI_Irecv(&dxm[0], 1, MPI_DOUBLE, comm_1d.west_rank, 2, comm_1d.mpi_comm, &request1[1]);
        } else {
            dxm[0] = dxm_f[0];
        }

        if (comm_1d.east_rank != MPI_PROC_NULL) {
            MPI_Isend(&dxm[nx],   1, MPI_DOUBLE, comm_1d.east_rank, 2, comm_1d.mpi_comm, &request2[0]);
            MPI_Irecv(&dxm[nx+1], 1, MPI_DOUBLE, comm_1d.east_rank, 1, comm_1d.mpi_comm, &request2[1]);
        } else {
            dxm[nx+1] = dxm_f[2*nx+1];
        }

        if (comm_1d.west_rank != MPI_PROC_NULL) {
            MPI_Waitall(2, request1, MPI_STATUSES_IGNORE);
        }
        if (comm_1d.east_rank != MPI_PROC_NULL) {
            MPI_Waitall(2, request2, MPI_STATUSES_IGNORE);
        }

        xg[0] = ox - 0.5 * dxm[0];
        for (i = 1; i <= nx; i++) {
            dxg[i] = 0.5 * (dxm[i] + dxm[i-1]);
            xg[i]  = xg[i-1] + dxg[i];
        }
        dxg[nx+1] = 0.5 * (dxm[nx+1] + dxm[nx]);
        xg[nx+1]  = xg[nx] + dxg[nx+1];
    } else {
        for (i = 0; i <= nx+1; i++) {
            dxm[i] = dxm_f[i];
            dxg[i] = dxg_f[i];
            xg[i]  = xg_f[i];
        }
    }
}

void multigrid_subdomain_aggregation_make_grid(double *dxm, double *dxg, double *xg, int nx, const double *dxm_f, const double *dxg_f, const double *xg_f, double ox, int lv_cur, int lv_coarsest, int lv_aggregation, int nx_lv_aggregation, cart_comm_1d comm_1d, char dir)
{
    int i;
    MPI_Request request1[2], request2[2];
    double *dxm_gl;

    if (lv_cur == lv_aggregation) {
        if (myrank == 0) {
            printf("[MG] level(aggregation) : %d, dir : %c\n", lv_cur, dir);
        }
    } else {
        if (myrank == 0) {
            printf("[MG] level : %d, dir : %c\n", lv_cur, dir);
        }
    } 

    for (i = 0; i <= nx+1; i++) {
        dxm[i] = 0.0;
        dxg[i] = 0.0;
        xg[i]  = 0.0;
    }

    if (lv_cur <= lv_coarsest) {
        if(lv_cur < lv_aggregation){
            for (i = 1; i <= nx; i++) {
                dxm[i] = dxm_f[2*i-1] + dxm_f[2*i];
            }

            if (comm_1d.west_rank != MPI_PROC_NULL) {
                MPI_Isend(&dxm[1], 1, MPI_DOUBLE, comm_1d.west_rank, 1, comm_1d.mpi_comm, &request1[0]);
                MPI_Irecv(&dxm[0], 1, MPI_DOUBLE, comm_1d.west_rank, 2, comm_1d.mpi_comm, &request1[1]);
            } else {
                dxm[0] = dxm_f[0];
            }

            if (comm_1d.east_rank != MPI_PROC_NULL) {
                MPI_Isend(&dxm[nx],   1, MPI_DOUBLE, comm_1d.east_rank, 2, comm_1d.mpi_comm, &request2[0]);
                MPI_Irecv(&dxm[nx+1], 1, MPI_DOUBLE, comm_1d.east_rank, 1, comm_1d.mpi_comm, &request2[1]);
            } else {
                dxm[nx+1] = dxm_f[2*nx+1];
            }

            if (comm_1d.west_rank != MPI_PROC_NULL) {
                MPI_Waitall(2, request1, MPI_STATUSES_IGNORE);
            }
            if (comm_1d.east_rank != MPI_PROC_NULL) {
                MPI_Waitall(2, request2, MPI_STATUSES_IGNORE);
            }

            xg[0] = ox - 0.5 * dxm[0];
            for (i = 1; i <= nx; i++) {
                dxg[i] = 0.5 * (dxm[i] + dxm[i-1]);
                xg[i]  = xg[i-1] + dxg[i];
            }
            dxg[nx+1] = 0.5 * (dxm[nx+1] + dxm[nx]);
            xg[nx+1]  = xg[nx] + dxg[nx+1];
        }
        else if(lv_cur == lv_aggregation){
            dxm_gl = (double*)malloc((nx_lv_aggregation+1) * sizeof(double));
            for (i = 1; i <= nx_lv_aggregation; i++) {
                dxm_gl[i] = dxm_f[2*i-1] + dxm_f[2*i];
            }
            MPI_Allgather(&dxm_gl[1], nx_lv_aggregation, MPI_DOUBLE, &dxm[1], nx_lv_aggregation, MPI_DOUBLE, comm_1d.mpi_comm);
            dxm[0]     = dxm_f[0];
            dxm[nx+1]  = dxm_f[2*nx_lv_aggregation + 1];
            MPI_Bcast(&dxm[0],    1, MPI_DOUBLE, 0, comm_1d.mpi_comm);
            MPI_Bcast(&dxm[nx+1], 1, MPI_DOUBLE, comm_1d.nprocs - 1, comm_1d.mpi_comm);

            xg[0] = ox - 0.5 * dxm[0];
            MPI_Bcast(&xg[0], 1, MPI_DOUBLE, 0, comm_1d.mpi_comm);
            for (i = 1; i <= nx; i++) {
                dxg[i] = 0.5 * (dxm[i] + dxm[i-1]);
                xg[i]  = xg[i-1] + dxg[i];
            }
            dxg[nx+1] = 0.5 * (dxm[nx+1] + dxm[nx]);
            xg[nx+1]  = xg[nx] + dxg[nx+1];
            free(dxm_gl);
        }
        else if(lv_cur > lv_aggregation){
            for (i = 1; i <= nx; i++) {
                dxm[i] = dxm_f[2*i-1] + dxm_f[2*i];
            }
            dxm[0]     = dxm_f[0];
            dxm[nx+1]  = dxm_f[2*nx + 1];
            
            xg[0] =  (xg_f[0] + 0.5 * dxm_f[0]) - 0.5 * dxm[0];
            for (i = 1; i <= nx; i++) {
                dxg[i] = 0.5 * (dxm[i] + dxm[i-1]);
                xg[i]  = xg[i-1] + dxg[i];
            }
            dxg[nx+1] = 0.5 * (dxm[nx+1] + dxm[nx]);
            xg[nx+1]  = xg[nx] + dxg[nx+1];
        }
    } 
    else {
        for (i = 0; i <= nx+1; i++) {
            dxm[i] = dxm_f[i];
            dxg[i] = dxg_f[i];
            xg[i]  = xg_f[i];
        }
    }
}

void multigrid_allocate_subdomain_variables() 
{
    int l;
    for (l = 1; l <= n_levels; l++) { 
        int nx = mg_sdm[l].nx;
        int ny = mg_sdm[l].ny;
        int nz = mg_sdm[l].nz;

        size_t total_size = (nx+2)*(ny+2)*(nz+2);
        mg_sdm[l].b = calloc(total_size, sizeof(double));
        if (!mg_sdm[l].b) { perror("calloc failed for b"); exit(1); }
        mg_sdm[l].r = calloc(total_size, sizeof(double));
        if (!mg_sdm[l].r) { perror("calloc failed for r"); exit(1); }
        mg_sdm[l].x = calloc(total_size, sizeof(double));
        if (!mg_sdm[l].x) { perror("calloc failed for x"); exit(1); }

    }

    for (l = 1; l <= n_levels; l++) {
        geometry_subdomain_ddt_create(&mg_sdm[l]);
    }
}

void multigrid_destroy() 
{
    int l;

    for (l = 1; l <= n_levels; l++) {  
        matrix_poisson_destroy(&mg_a_poisson[l]);
        geometry_subdomain_destroy(&mg_sdm[l]);
        geometry_subdomain_ddt_destroy(&mg_sdm[l]);
    }

    free(mg_sdm);
    mg_sdm = NULL;

    free(mg_a_poisson);
    mg_a_poisson = NULL;
}

void multigrid_restriction(double *val_c, double *val_f, subdomain *dm_c, subdomain *dm_f, int level) 
{
    int i, j, k;
    int i_c, j_c, k_c;
    int iz_f, jz_f, kz_f;
    int ip_f, jp_f, kp_f;
    int i_gl_c, j_gl_c, k_gl_c;
    int i_stride_f, i_offset_f;
    int j_stride_f, j_offset_f;
    int k_stride_f, k_offset_f;
    int nx_c, ny_c, nz_c;
    int nx_f, ny_f, nz_f;
    double vol_f[2][2][2], vol_c;
    double *dxf = (double*)calloc(dm_f->nx + 2, sizeof(double));
    double *dyf = (double*)calloc(dm_f->ny + 2, sizeof(double));
    double *dzf = (double*)calloc(dm_f->nz + 2, sizeof(double));

    nx_f = dm_f->nx;
    ny_f = dm_f->ny;
    nz_f = dm_f->nz;

    switch (aggregation_type)
    {
        case 0:
            nx_c = dm_c->nx; ny_c = dm_c->ny; nz_c = dm_c->nz;
            i_gl_c = j_gl_c = k_gl_c = 0;
            break;
        case 1:
            if (level == lv_aggregation-1) {
                nx_c = dm_c->nx / comm_1d_x.nprocs;
                ny_c = dm_c->ny / comm_1d_y.nprocs;
                nz_c = dm_c->nz / comm_1d_z.nprocs;
                i_gl_c = nx_c * comm_1d_x.myrank;
                j_gl_c = ny_c * comm_1d_y.myrank;
                k_gl_c = nz_c * comm_1d_z.myrank;
            } else {
                nx_c = dm_c->nx; ny_c = dm_c->ny; nz_c = dm_c->nz;
                i_gl_c = j_gl_c = k_gl_c = 0;
            }
            break;
        case 2:
            nx_c = (level == lv_aggregation_x-1) ? dm_c->nx / comm_1d_x.nprocs : dm_c->nx;
            i_gl_c = (level == lv_aggregation_x-1) ? nx_c * comm_1d_x.myrank : 0;
            ny_c = (level == lv_aggregation_y-1) ? dm_c->ny / comm_1d_y.nprocs : dm_c->ny;
            j_gl_c = (level == lv_aggregation_y-1) ? ny_c * comm_1d_y.myrank : 0;
            nz_c = (level == lv_aggregation_z-1) ? dm_c->nz / comm_1d_z.nprocs : dm_c->nz;
            k_gl_c = (level == lv_aggregation_z-1) ? nz_c * comm_1d_z.myrank : 0;
            break;
        default:
            if (myrank == 0) printf("[Error] Aggregation method should be 0, 1, or 2: %d\n", aggregation_type);
            exit(1);
    }

    if (level >= lv_gdm_coarsest_x) {
        for (i = 0; i <= dm_f->nx+1; i++) dxf[i] = dm_f->dxm[i];
        i_stride_f = 1; 
        i_offset_f = 0;
    } else {
        for (i = 1; i <= dm_f->nx; i++) {
            i_c = (i-1)/2 + 1;
            dxf[i] = dm_c->dxm[i_c] - dm_f->dxm[i];
        }
        i_stride_f = 2;
        i_offset_f = 1;
    }

    if (level >= lv_gdm_coarsest_y) {
        for (j = 0; j <= dm_f->ny+1; j++) dyf[j] = dm_f->dym[j];
        j_stride_f = 1; 
        j_offset_f = 0;
    } else {
        for (j = 1; j <= dm_f->ny; j++) {
            j_c = (j-1)/2 + 1;
            dyf[j] = dm_c->dym[j_c] - dm_f->dym[j];
        }
        j_stride_f = 2; 
        j_offset_f = 1;
    }

    if (level >= lv_gdm_coarsest_z) {
        for (k = 0; k <= dm_f->nz+1; k++) dzf[k] = dm_f->dzm[k];
        k_stride_f = 1; 
        k_offset_f = 0;
    } else {
        for (k = 1; k <= dm_f->nz; k++) {
            k_c = (k-1)/2 + 1;
            dzf[k] = dm_c->dzm[k_c] - dm_f->dzm[k];
        }
        k_stride_f = 2; 
        k_offset_f = 1;
    }

    
    for (i = 1; i <= nx_c; i++) 
    {
        i_c = i + i_gl_c;
        ip_f = i * i_stride_f; 
        iz_f = ip_f - i_offset_f;
        for (j = 1; j <= ny_c; j++) 
        {
            j_c = j + j_gl_c;
            jp_f = j * j_stride_f; 
            jz_f = jp_f - j_offset_f;
            for (k = 1; k <= nz_c; k++)
            {
                k_c = k + k_gl_c;
                kp_f = k * k_stride_f; 
                kz_f = kp_f - k_offset_f;
        
                vol_f[0][0][0] = dxf[iz_f] * dyf[jz_f] * dzf[kz_f];
                vol_f[1][0][0] = dxf[ip_f] * dyf[jz_f] * dzf[kz_f];
                vol_f[0][1][0] = dxf[iz_f] * dyf[jp_f] * dzf[kz_f];
                vol_f[1][1][0] = dxf[ip_f] * dyf[jp_f] * dzf[kz_f];
                vol_f[0][0][1] = dxf[iz_f] * dyf[jz_f] * dzf[kp_f];
                vol_f[1][0][1] = dxf[ip_f] * dyf[jz_f] * dzf[kp_f];
                vol_f[0][1][1] = dxf[iz_f] * dyf[jp_f] * dzf[kp_f];
                vol_f[1][1][1] = dxf[ip_f] * dyf[jp_f] * dzf[kp_f];

                vol_c = 0.0;
                for (int ii=0; ii<2; ii++)
                {
                    for (int jj=0; jj<2; jj++)
                    {
                        for (int kk=0; kk<2; kk++)
                        {
                            vol_c += vol_f[ii][jj][kk];
                        }
                            
                    }
                        
                }

                val_c[IDX(i_c,j_c,k_c,dm_c->nx,dm_c->ny,dm_c->nz)] = ( vol_f[0][0][0] * val_f[IDX(ip_f,jp_f,kp_f,nx_f,ny_f,nz_f)] 
                                                         + vol_f[1][0][0] * val_f[IDX(iz_f,jp_f,kp_f,nx_f,ny_f,nz_f)] 
                                                         + vol_f[0][1][0] * val_f[IDX(ip_f,jz_f,kp_f,nx_f,ny_f,nz_f)] 
                                                         + vol_f[1][1][0] * val_f[IDX(iz_f,jz_f,kp_f,nx_f,ny_f,nz_f)] 
                                                         + vol_f[0][0][1] * val_f[IDX(ip_f,jp_f,kz_f,nx_f,ny_f,nz_f)] 
                                                         + vol_f[1][0][1] * val_f[IDX(iz_f,jp_f,kz_f,nx_f,ny_f,nz_f)] 
                                                         + vol_f[0][1][1] * val_f[IDX(ip_f,jz_f,kz_f,nx_f,ny_f,nz_f)] 
                                                         + vol_f[1][1][1] * val_f[IDX(iz_f,jz_f,kz_f,nx_f,ny_f,nz_f)] ) / vol_c;
            }
        }
    }
    free(dxf); free(dyf); free(dzf);
    dxf = dyf = dzf = NULL;

#ifdef DEBUG_RESTRICTION
    multigrid_common_print_restriction_info(val_c, val_f, vol_c, vol_f, level, nx_c, ny_c, nz_c, i_gl_c, j_gl_c, k_gl_c, &
                                                    i_stride_f, j_stride_f, k_stride_f, i_offset_f, j_offset_f, k_offset_f)
#endif
                

}





void multigrid_prolongation_linear_on_nonuniform_grid(double *val_f, const double *val_c, subdomain *dm_f, subdomain *dm_c, int level)
{
    int nx_c, ny_c, nz_c;
    int nx_f, ny_f, nz_f;
    int i, j, k;
    int iz_f, jz_f, kz_f;
    int ip_f, jp_f, kp_f;
    int i_stride_f, i_offset_f;
    int j_stride_f, j_offset_f;
    int k_stride_f, k_offset_f;
    int i_gl_c, j_gl_c, k_gl_c;
    int im_c, jm_c, km_c;
    int iz_c, jz_c, kz_c;
    int ip_c, jp_c, kp_c;

    nx_f = dm_f->nx;
    ny_f = dm_f->ny;
    nz_f = dm_f->nz;

    double *dxp = (double*)calloc(dm_f->nx+2, sizeof(double));
    double *dxn = (double*)calloc(dm_f->nx+2, sizeof(double));
    double *dyp = (double*)calloc(dm_f->ny+2, sizeof(double));
    double *dyn = (double*)calloc(dm_f->ny+2, sizeof(double));
    double *dzp = (double*)calloc(dm_f->nz+2, sizeof(double));
    double *dzn = (double*)calloc(dm_f->nz+2, sizeof(double));



    switch(aggregation_type) {
        case 0:
            nx_c = dm_c->nx; ny_c = dm_c->ny; nz_c = dm_c->nz;
            i_gl_c = j_gl_c = k_gl_c = 0;
            break;
        case 1:
            if(level == lv_aggregation-1) {
                nx_c = dm_c->nx / comm_1d_x.nprocs;
                ny_c = dm_c->ny / comm_1d_y.nprocs;
                nz_c = dm_c->nz / comm_1d_z.nprocs;
                i_gl_c = nx_c * comm_1d_x.myrank;
                j_gl_c = ny_c * comm_1d_y.myrank;
                k_gl_c = nz_c * comm_1d_z.myrank;
            } else {
                nx_c = dm_c->nx; ny_c = dm_c->ny; nz_c = dm_c->nz;
                i_gl_c = j_gl_c = k_gl_c = 0;
            }
            break;
        case 2:
            nx_c = (level == lv_aggregation_x-1) ? dm_c->nx / comm_1d_x.nprocs : dm_c->nx;
            i_gl_c = (level == lv_aggregation_x-1) ? nx_c * comm_1d_x.myrank : 0;
            ny_c = (level == lv_aggregation_y-1) ? dm_c->ny / comm_1d_y.nprocs : dm_c->ny;
            j_gl_c = (level == lv_aggregation_y-1) ? ny_c * comm_1d_y.myrank : 0;
            nz_c = (level == lv_aggregation_z-1) ? dm_c->nz / comm_1d_z.nprocs : dm_c->nz;
            k_gl_c = (level == lv_aggregation_z-1) ? nz_c * comm_1d_z.myrank : 0;
            break;
        default:
            if(myrank==0) printf("[Error] Aggregation method should be 0,1,2: %d\n", aggregation_type);
            exit(1);
    }

   

    if(level >= lv_gdm_coarsest_x) {
        i_stride_f = 1; i_offset_f = 0;
        for(i=0;i<=dm_c->nx+1;i++){ dxp[i]=dxn[i]=dm_c->dxm[i]; }
    } else {
        i_stride_f = 2; i_offset_f = 1;
        for(i=1;i<=nx_c;i++){
            ip_f = 2*i;
            dxp[ip_f-1] = 0.5*dm_c->dxm[i-1+i_gl_c] + 0.5*dm_f->dxm[ip_f-1];
            dxp[ip_f]   = 0.5*dm_c->dxm[i+i_gl_c]   - 0.5*dm_f->dxm[ip_f-1];
            dxn[ip_f-1] = 0.5*dm_c->dxm[i+i_gl_c]   - 0.5*dm_f->dxm[ip_f];
            dxn[ip_f]   = 0.5*dm_c->dxm[i+1+i_gl_c] + 0.5*dm_f->dxm[ip_f];
        }
    }

    if(level >= lv_gdm_coarsest_y) {
        j_stride_f = 1; j_offset_f = 0;
        for(j=0;j<=dm_c->ny+1;j++){ dyp[j]=dyn[j]=dm_c->dym[j]; }
    } else {
        j_stride_f = 2; j_offset_f = 1;
        for(j=1;j<=ny_c;j++){
            jp_f = 2*j;
            dyp[jp_f-1] = 0.5*dm_c->dym[j-1+j_gl_c] + 0.5*dm_f->dym[jp_f-1];
            dyp[jp_f]   = 0.5*dm_c->dym[j+j_gl_c]   - 0.5*dm_f->dym[jp_f-1];
            dyn[jp_f-1] = 0.5*dm_c->dym[j+j_gl_c]   - 0.5*dm_f->dym[jp_f];
            dyn[jp_f]   = 0.5*dm_c->dym[j+1+j_gl_c] + 0.5*dm_f->dym[jp_f];
        }
    }

    if(level >= lv_gdm_coarsest_z) {
        k_stride_f = 1; k_offset_f = 0;
        for(k=0;k<=dm_c->nz+1;k++){ dzp[k]=dzn[k]=dm_c->dzm[k]; }
    } else {
        k_stride_f = 2; k_offset_f = 1;
        for(k=1;k<=nz_c;k++){
            kp_f = 2*k;
            dzp[kp_f-1] = 0.5*dm_c->dzm[k-1+k_gl_c] + 0.5*dm_f->dzm[kp_f-1];
            dzp[kp_f]   = 0.5*dm_c->dzm[k+k_gl_c]   - 0.5*dm_f->dzm[kp_f-1];
            dzn[kp_f-1] = 0.5*dm_c->dzm[k+k_gl_c]   - 0.5*dm_f->dzm[kp_f];
            dzn[kp_f]   = 0.5*dm_c->dzm[k+1+k_gl_c] + 0.5*dm_f->dzm[kp_f];
        }
    }
    
    
    for (i=1;i<=nx_c;i++)
    {
        ip_f = i * i_stride_f;
        iz_f = ip_f - i_offset_f;
        iz_c = i + i_gl_c;
        im_c = iz_c - i_offset_f;
        ip_c = iz_c + i_offset_f;

        for(j=1;j<=ny_c;j++)
        {
            jp_f = j * j_stride_f;
            jz_f = jp_f - j_offset_f;
            jz_c = j + j_gl_c;
            jm_c = jz_c - j_offset_f;
            jp_c = jz_c + j_offset_f;

            for (k=1;k<=nz_c;k++)
            {
                kp_f = k * k_stride_f;
                kz_f = kp_f - k_offset_f;
                kz_c = k + k_gl_c;
                km_c = kz_c - k_offset_f;
                kp_c = kz_c + k_offset_f;

                

                val_f[IDX(iz_f,jz_f,kz_f,nx_f,ny_f,nz_f)] = (dxp[iz_f]*dyp[jz_f]*dzp[kz_f]*val_c[IDX(iz_c,jz_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxp[ip_f]*dyp[jz_f]*dzp[kz_f]*val_c[IDX(im_c,jz_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxp[iz_f]*dyp[jp_f]*dzp[kz_f]*val_c[IDX(iz_c,jm_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxp[iz_f]*dyp[jz_f]*dzp[kp_f]*val_c[IDX(iz_c,jz_c,km_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxp[ip_f]*dyp[jp_f]*dzp[kz_f]*val_c[IDX(im_c,jm_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxp[ip_f]*dyp[jz_f]*dzp[kp_f]*val_c[IDX(im_c,jz_c,km_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxp[iz_f]*dyp[jp_f]*dzp[kp_f]*val_c[IDX(iz_c,jm_c,km_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxp[ip_f]*dyp[jp_f]*dzp[kp_f]*val_c[IDX(im_c,jm_c,km_c,dm_c->nx,dm_c->ny,dm_c->nz)] )
                                                        / ( (dxp[iz_f]+dxp[ip_f])*(dyp[jz_f]+dyp[jp_f])*(dzp[kz_f]+dzp[kp_f]) );
               
                val_f[IDX(ip_f,jz_f,kz_f,nx_f,ny_f,nz_f)] = (dxn[ip_f]*dyp[jz_f]*dzp[kz_f]*val_c[IDX(iz_c,jz_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxn[iz_f]*dyp[jz_f]*dzp[kz_f]*val_c[IDX(ip_c,jz_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxn[ip_f]*dyp[jp_f]*dzp[kz_f]*val_c[IDX(iz_c,jm_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxn[ip_f]*dyp[jz_f]*dzp[kp_f]*val_c[IDX(iz_c,jz_c,km_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxn[iz_f]*dyp[jp_f]*dzp[kz_f]*val_c[IDX(ip_c,jm_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxn[iz_f]*dyp[jz_f]*dzp[kp_f]*val_c[IDX(ip_c,jz_c,km_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxn[ip_f]*dyp[jp_f]*dzp[kp_f]*val_c[IDX(iz_c,jm_c,km_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxn[iz_f]*dyp[jp_f]*dzp[kp_f]*val_c[IDX(ip_c,jm_c,km_c,dm_c->nx,dm_c->ny,dm_c->nz)] )
                                                        / ( (dxn[iz_f]+dxn[ip_f])*(dyp[jz_f]+dyp[jp_f])*(dzp[kz_f]+dzp[kp_f]) );
              
                val_f[IDX(iz_f,jp_f,kz_f,nx_f,ny_f,nz_f)] = (dxp[iz_f]*dyn[jp_f]*dzp[kz_f]*val_c[IDX(iz_c,jz_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxp[ip_f]*dyn[jp_f]*dzp[kz_f]*val_c[IDX(im_c,jz_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxp[iz_f]*dyn[jz_f]*dzp[kz_f]*val_c[IDX(iz_c,jp_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxp[iz_f]*dyn[jp_f]*dzp[kp_f]*val_c[IDX(iz_c,jz_c,km_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxp[ip_f]*dyn[jz_f]*dzp[kz_f]*val_c[IDX(im_c,jp_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxp[ip_f]*dyn[jp_f]*dzp[kp_f]*val_c[IDX(im_c,jz_c,km_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxp[iz_f]*dyn[jz_f]*dzp[kp_f]*val_c[IDX(iz_c,jp_c,km_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxp[ip_f]*dyn[jz_f]*dzp[kp_f]*val_c[IDX(im_c,jp_c,km_c,dm_c->nx,dm_c->ny,dm_c->nz)] )
                                                        / ( (dxp[iz_f]+dxp[ip_f])*(dyn[jz_f]+dyn[jp_f])*(dzp[kz_f]+dzp[kp_f]) );
             
                val_f[IDX(ip_f,jp_f,kz_f,nx_f,ny_f,nz_f)] = (dxn[ip_f]*dyn[jp_f]*dzp[kz_f]*val_c[IDX(iz_c,jz_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxn[iz_f]*dyn[jp_f]*dzp[kz_f]*val_c[IDX(ip_c,jz_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxn[ip_f]*dyn[jz_f]*dzp[kz_f]*val_c[IDX(iz_c,jp_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxn[ip_f]*dyn[jp_f]*dzp[kp_f]*val_c[IDX(iz_c,jz_c,km_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxn[iz_f]*dyn[jz_f]*dzp[kz_f]*val_c[IDX(ip_c,jp_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxn[iz_f]*dyn[jp_f]*dzp[kp_f]*val_c[IDX(ip_c,jz_c,km_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxn[ip_f]*dyn[jz_f]*dzp[kp_f]*val_c[IDX(iz_c,jp_c,km_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxn[iz_f]*dyn[jz_f]*dzp[kp_f]*val_c[IDX(ip_c,jp_c,km_c,dm_c->nx,dm_c->ny,dm_c->nz)] )
                                                        / ( (dxn[iz_f]+dxn[ip_f])*(dyn[jz_f]+dyn[jp_f])*(dzp[kz_f]+dzp[kp_f]) );
                val_f[IDX(iz_f,jz_f,kp_f,nx_f,ny_f,nz_f)] = (dxp[iz_f]*dyp[jz_f]*dzn[kp_f]*val_c[IDX(iz_c,jz_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxp[ip_f]*dyp[jz_f]*dzn[kp_f]*val_c[IDX(im_c,jz_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxp[iz_f]*dyp[jp_f]*dzn[kp_f]*val_c[IDX(iz_c,jm_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxp[iz_f]*dyp[jz_f]*dzn[kz_f]*val_c[IDX(iz_c,jz_c,kp_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxp[ip_f]*dyp[jp_f]*dzn[kp_f]*val_c[IDX(im_c,jm_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxp[ip_f]*dyp[jz_f]*dzn[kz_f]*val_c[IDX(im_c,jz_c,kp_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxp[iz_f]*dyp[jp_f]*dzn[kz_f]*val_c[IDX(iz_c,jm_c,kp_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxp[ip_f]*dyp[jp_f]*dzn[kz_f]*val_c[IDX(im_c,jm_c,kp_c,dm_c->nx,dm_c->ny,dm_c->nz)] )
                                                        / ( (dxp[iz_f]+dxp[ip_f])*(dyp[jz_f]+dyp[jp_f])*(dzn[kz_f]+dzn[kp_f]) );
                val_f[IDX(ip_f,jz_f,kp_f,nx_f,ny_f,nz_f)] = (dxn[ip_f]*dyp[jz_f]*dzn[kp_f]*val_c[IDX(iz_c,jz_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxn[iz_f]*dyp[jz_f]*dzn[kp_f]*val_c[IDX(ip_c,jz_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxn[ip_f]*dyp[jp_f]*dzn[kp_f]*val_c[IDX(iz_c,jm_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxn[ip_f]*dyp[jz_f]*dzn[kz_f]*val_c[IDX(iz_c,jz_c,kp_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxn[iz_f]*dyp[jp_f]*dzn[kp_f]*val_c[IDX(ip_c,jm_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxn[iz_f]*dyp[jz_f]*dzn[kz_f]*val_c[IDX(ip_c,jz_c,kp_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxn[ip_f]*dyp[jp_f]*dzn[kz_f]*val_c[IDX(iz_c,jm_c,kp_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxn[iz_f]*dyp[jp_f]*dzn[kz_f]*val_c[IDX(ip_c,jm_c,kp_c,dm_c->nx,dm_c->ny,dm_c->nz)] )
                                                        / ( (dxn[iz_f]+dxn[ip_f])*(dyp[jz_f]+dyp[jp_f])*(dzn[kz_f]+dzn[kp_f]) );
                val_f[IDX(iz_f,jp_f,kp_f,nx_f,ny_f,nz_f)] = (dxp[iz_f]*dyn[jp_f]*dzn[kp_f]*val_c[IDX(iz_c,jz_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxp[ip_f]*dyn[jp_f]*dzn[kp_f]*val_c[IDX(im_c,jz_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxp[iz_f]*dyn[jz_f]*dzn[kp_f]*val_c[IDX(iz_c,jp_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxp[iz_f]*dyn[jp_f]*dzn[kz_f]*val_c[IDX(iz_c,jz_c,kp_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxp[ip_f]*dyn[jz_f]*dzn[kp_f]*val_c[IDX(im_c,jp_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxp[ip_f]*dyn[jp_f]*dzn[kz_f]*val_c[IDX(im_c,jz_c,kp_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxp[iz_f]*dyn[jz_f]*dzn[kz_f]*val_c[IDX(iz_c,jp_c,kp_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxp[ip_f]*dyn[jz_f]*dzn[kz_f]*val_c[IDX(im_c,jp_c,kp_c,dm_c->nx,dm_c->ny,dm_c->nz)] )
                                                        / ( (dxp[iz_f]+dxp[ip_f])*(dyn[jz_f]+dyn[jp_f])*(dzn[kz_f]+dzn[kp_f]) );
                val_f[IDX(ip_f,jp_f,kp_f,nx_f,ny_f,nz_f)] = (dxn[ip_f]*dyn[jp_f]*dzn[kp_f]*val_c[IDX(iz_c,jz_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxn[iz_f]*dyn[jp_f]*dzn[kp_f]*val_c[IDX(ip_c,jz_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxn[ip_f]*dyn[jz_f]*dzn[kp_f]*val_c[IDX(iz_c,jp_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxn[ip_f]*dyn[jp_f]*dzn[kz_f]*val_c[IDX(iz_c,jz_c,kp_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxn[iz_f]*dyn[jz_f]*dzn[kp_f]*val_c[IDX(ip_c,jp_c,kz_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxn[iz_f]*dyn[jp_f]*dzn[kz_f]*val_c[IDX(ip_c,jz_c,kp_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxn[ip_f]*dyn[jz_f]*dzn[kz_f]*val_c[IDX(iz_c,jp_c,kp_c,dm_c->nx,dm_c->ny,dm_c->nz)] 
                                                            +dxn[iz_f]*dyn[jz_f]*dzn[kz_f]*val_c[IDX(ip_c,jp_c,kp_c,dm_c->nx,dm_c->ny,dm_c->nz)] )
                                                        / ( (dxn[iz_f]+dxn[ip_f])*(dyn[jz_f]+dyn[jp_f])*(dzn[kz_f]+dzn[kp_f]) );

            
                    
            }
        }
    }
    
    free(dxp); free(dxn);
    free(dyp); free(dyn);
    free(dzp); free(dzn);
    dxp = dxn = NULL;
    dyp = dyn = NULL;
    dzp = dzn = NULL;

#ifdef DEBUG_PROLONGATION
        multigrid_common_print_prolongation_info(val_c, val_f, level, dm_f%nx, dm_f%ny, dm_f%nz, i_gl_c, j_gl_c, k_gl_c, &
                                                    i_stride_f, j_stride_f, k_stride_f, i_offset_f, j_offset_f, k_offset_f)
#endif

}

void multigrid_residual(double *rsd, matrix_poisson *a_poisson, double *x, double *rhs, subdomain *dm, int is_aggregated[3])
{
    int i,j,k;
    int nx = dm->nx;
    int ny = dm->ny;
    int nz = dm->nz;

    #pragma omp parallel for collapse(3)
    for(i=0;i<=nx+1;i++){
        for(j=0;j<=ny+1;j++){
            for(k=0;k<=nz+1;k++){
                rsd[IDX(i,j,k,nx,ny,nz)] = 0.0;
            }
        }
    }

    mv_mul_poisson_matrix(rsd, a_poisson, x, dm, is_aggregated);

    #pragma omp parallel for collapse(3)
    for(i=1;i<=nx;i++){
        for(j=1;j<=ny;j++){
            for(k=1;k<=nz;k++){
                rsd[IDX(i,j,k,nx,ny,nz)] = rhs[IDX(i,j,k,nx,ny,nz)] - rsd[IDX(i,j,k,nx,ny,nz)];
            }
        }
    }
}

void multigrid_solve_coarset_level(double *x, matrix_poisson *a_poisson, double *rhs, subdomain *dm, int maxiteration, double tolerance, double omega, int is_aggregated[3])
{
    int i, j, k;
    int nx = dm->nx;
    int ny = dm->ny;
    int nz = dm->nz;

    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int size = (nx+2)*(ny+2)*(nz+2);
    for (i=0; i<size; i++) x[i] = 0.0;

    if (myrank == 0) {
        printf("[MG] Obtain the solution on a grid in the coarsest level\n");
    }
  

    if ((nx*ny*nz == 1) && (is_aggregated[0] && is_aggregated[1] && is_aggregated[2])) {
        for (i=1; i<=nx; i++) {
            for (j=1; j<=ny; j++) {
                for (k=1; k<=nz; k++) {
                    size_t idx = IDX(i,j,k,nx,ny,nz);
                    x[idx] = rhs[idx] / *COEFF(a_poisson, 0, i, j, k);
                }
            }
        }
    } else {
        rbgs_solver_poisson_matrix(x, a_poisson, rhs, dm, maxiteration, tolerance, omega, is_aggregated);
    }
}

void multigrid_solve_vcycle(double *sol, matrix_poisson *a_poisson, double *rhs, subdomain *sdm, int maxiteration, double tolerance, double omega_sor)
{
    int myrank = -1;
    int mpi_inited = 0;
    MPI_Initialized(&mpi_inited);
    if (mpi_inited) MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    switch (aggregation_type) {
        case 0:
            multigrid_no_aggregation_vcycle_solver(sol, a_poisson, rhs, sdm, maxiteration, tolerance, omega_sor);
            break;
        case 1:
            multigrid_single_aggregation_vcycle_solver(sol, a_poisson, rhs, sdm, maxiteration, tolerance, omega_sor);
            break;
        case 2:
            multigrid_adaptive_aggregation_vcycle_solver(sol, a_poisson, rhs, sdm, maxiteration, tolerance, omega_sor);
            break;
        default: {
            if (myrank == 0) {
                fprintf(stderr, "[Error] Aggregation method should be 0, 1, or 2. Got %d\n", aggregation_type);
            }
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
    }
}

void multigrid_no_aggregation_vcycle_solver(double *sol, matrix_poisson *a_poisson, double *rhs, subdomain *sdm, int maxiteration, double tolerance, double omega_sor)
{
    int l, cyc, nx, ny, nz;
    double rsd_val, res0tol;

    nx = sdm->nx;
    ny = sdm->ny;
    nz = sdm->nz;
    size_t size = (nx+2)*(ny+2)*(nz+2);
    double *rsd = (double*) calloc(size, sizeof(double));

    multigrid_residual(rsd, a_poisson, sol, rhs, sdm, sdm->is_aggregated);
    vv_dot_3d_matrix(&res0tol, rsd, rsd, nx, ny, nz, sdm->is_aggregated);

    for(cyc = 1; cyc <= n_vcycles; cyc++)
    {
        rbgs_iterator_poisson_matrix(sol, a_poisson, rhs, sdm, maxiteration, omega_sor, sdm->is_aggregated);
        multigrid_residual(rsd, a_poisson, sol, rhs, sdm, sdm->is_aggregated);
        multigrid_restriction(mg_sdm[1].b, rsd, &mg_sdm[1], sdm, 0);
        if(myrank==0) printf("[MG] Restriction from level 0 to level 1\n");

        for(l=1; l<=n_levels-1; l++)
        {
            size_t level_size = (mg_sdm[l].nx + 2) * (mg_sdm[l].ny + 2) * (mg_sdm[l].nz + 2);
            memset(mg_sdm[l].x, 0, level_size * sizeof(double));
            rbgs_iterator_poisson_matrix(mg_sdm[l].x, &mg_a_poisson[l], mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
            multigrid_residual(mg_sdm[l].r, &mg_a_poisson[l], mg_sdm[l].x, mg_sdm[l].b, &mg_sdm[l], mg_sdm[l].is_aggregated);
            multigrid_restriction(mg_sdm[l+1].b, mg_sdm[l].r, &mg_sdm[l+1], &mg_sdm[l], l);
            if(myrank==0) printf("[MG] Restriction from level %d to level %d\n", l, l+1);
        }

        multigrid_solve_coarset_level(mg_sdm[n_levels].x, &mg_a_poisson[n_levels], mg_sdm[n_levels].b, &mg_sdm[n_levels], 1000, tolerance, omega_sor, mg_sdm[n_levels].is_aggregated);

        multigrid_residual(mg_sdm[n_levels].r, &mg_a_poisson[n_levels], mg_sdm[n_levels].x, mg_sdm[n_levels].b, &mg_sdm[n_levels], mg_sdm[n_levels].is_aggregated);
        if (myrank == 0) {
            printf("[MG] Solution in the coarsest level : x(1,1,1) = %18.10e, residue = %18.10e\n",
            mg_sdm[n_levels].x[IDX(1,1,1, mg_sdm[n_levels].nx, mg_sdm[n_levels].ny, mg_sdm[n_levels].nz)],
            mg_sdm[n_levels].r[IDX(1,1,1, mg_sdm[n_levels].nx, mg_sdm[n_levels].ny, mg_sdm[n_levels].nz)]);
        }
    #ifdef DEBUG_COARSEST
        multigrid_common_print_coarsest_level_solution(cyc, &mg_sdm[n_levels]);
    #endif

        for(l = n_levels-1; l >= 1; l--)
        {
            geometry_halocell_update_selectively(mg_sdm[l+1].x, &mg_sdm[l+1], mg_sdm[l+1].is_aggregated);
            multigrid_prolongation_linear_on_nonuniform_grid(mg_sdm[l].r, mg_sdm[l+1].x, &mg_sdm[l], &mg_sdm[l+1], l);
            
            size_t level_size = (mg_sdm[l].nx + 2) * (mg_sdm[l].ny + 2) * (mg_sdm[l].nz + 2);
            for(size_t idx=0; idx<level_size; idx++) mg_sdm[l].x[idx] = mg_sdm[l].x[idx] + mg_sdm[l].r[idx];

            rbgs_iterator_poisson_matrix(mg_sdm[l].x, &mg_a_poisson[l], mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
        }

        geometry_halocell_update_selectively(mg_sdm[1].x, &mg_sdm[1], mg_sdm[1].is_aggregated);
        multigrid_prolongation_linear_on_nonuniform_grid(rsd, mg_sdm[1].x, sdm, &mg_sdm[1], 0);
        for(size_t idx=0; idx<size; idx++) sol[idx] = sol[idx] + rsd[idx];
        rbgs_iterator_poisson_matrix(sol, a_poisson, rhs, sdm, maxiteration, omega_sor, sdm->is_aggregated);
        multigrid_residual(rsd, a_poisson, sol, rhs, sdm, sdm->is_aggregated);
        vv_dot_3d_matrix(&rsd_val, rsd, rsd, nx, ny, nz, sdm->is_aggregated);

        if(myrank==0)
            printf("[MG] cycle = %d, Error = %e, Initial = %e, Relative = %e\n", cyc, sqrt(rsd_val), sqrt(res0tol), sqrt(rsd_val/res0tol));

        if(sqrt(rsd_val/res0tol) < tolerance) break;
    }

    if(myrank==0) printf("[MG] Total %d V-cycles end\n", cyc);

    free(rsd);
    rsd = NULL;

}

void multigrid_single_aggregation_vcycle_solver(double *sol, matrix_poisson *a_poisson, double *rhs, subdomain *sdm, int maxiteration, double tolerance, double omega_sor)
{
    int l, i, cyc;
    double rsd_val, res0tol;

    MPI_Datatype  ddtype_temp1, ddtype_send, ddtype_gatherv;
    MPI_Datatype ddtype_scatterv, ddtype_recv;
    int sizes[3], subsizes[3], starts[3];
    int r8size;
    int nx_aggr, ny_aggr, nz_aggr;
    int nx, ny, nz;
    int ix, iy, iz;
    MPI_Aint extent, lb;
    int *cnt_gatherv, *disps_gatherv, *cnt_scatterv, *disps_scatterv;
    int (*cart_coord)[3];
    cart_coord = malloc(nprocs * sizeof(*cart_coord));
    int ierr;
    for (i = 0; i < nprocs; i++) {
        ierr = MPI_Cart_coords(mpi_world_cart, i, 3, cart_coord[i]);
        if (ierr != MPI_SUCCESS) {
            fprintf(stderr, "MPI_Cart_coords error\n");
            MPI_Abort(MPI_COMM_WORLD, ierr);
        }
    }
    cnt_gatherv = malloc((nprocs ) * sizeof(int));
    disps_gatherv = malloc((nprocs ) * sizeof(int));

    nx_aggr = mg_sdm[lv_aggregation].nx;
    ny_aggr = mg_sdm[lv_aggregation].ny;
    nz_aggr = mg_sdm[lv_aggregation].nz;

    nx = nx_aggr / comm_1d_x.nprocs;
    ny = ny_aggr / comm_1d_y.nprocs;
    nz = nz_aggr / comm_1d_z.nprocs;

    ix = nx * comm_1d_x.myrank;
    iy = ny * comm_1d_y.myrank;
    iz = nz * comm_1d_z.myrank;

    // For MPI_Gatherv
    sizes[0] = nx_aggr + 2;
    sizes[1] = ny_aggr + 2;
    sizes[2] = nz_aggr + 2;
    subsizes[0] = nx;
    subsizes[1] = ny;
    subsizes[2] = nz;
    starts[0] = ix+1;
    starts[1] = iy+1;
    starts[2] = iz+1;
    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &ddtype_temp1);
    MPI_Type_size(MPI_DOUBLE, &r8size);
    lb = 0;
    extent = (MPI_Aint) r8size;
    MPI_Type_create_resized(ddtype_temp1, lb, extent, &ddtype_send);
    MPI_Type_commit(&ddtype_send);

    
    starts[0] = 1;
    starts[1] = 1;
    starts[2] = 1;
    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &ddtype_temp1);
    MPI_Type_size(MPI_DOUBLE, &r8size);
    lb = 0;
    extent = (MPI_Aint) r8size;
    MPI_Type_create_resized(ddtype_temp1, lb, extent, &ddtype_gatherv);
    MPI_Type_commit(&ddtype_gatherv);


    for (i = 0; i < nprocs; i++) {
        cnt_gatherv[i] = 1;
        disps_gatherv[i] = nz * (cart_coord[i][2])
                         + ny * (cart_coord[i][1]) * (nz_aggr + 2)       
                         + nx * (cart_coord[i][0]) * (nz_aggr + 2) * (ny_aggr + 2);    
    }

    // For MPI_Scatterv
    cnt_scatterv = malloc((nprocs ) * sizeof(int));
    disps_scatterv = malloc((nprocs ) * sizeof(int));

    nx_aggr = mg_sdm[lv_aggregation].nx;
    ny_aggr = mg_sdm[lv_aggregation].ny;
    nz_aggr = mg_sdm[lv_aggregation].nz;

    nx = nx_aggr / comm_1d_x.nprocs;
    ny = ny_aggr / comm_1d_y.nprocs;
    nz = nz_aggr / comm_1d_z.nprocs;

    ix = nx * comm_1d_x.myrank;
    iy = ny * comm_1d_y.myrank;
    iz = nz * comm_1d_z.myrank;

    sizes[0] = nx_aggr + 2;
    sizes[1] = ny_aggr + 2;
    sizes[2] = nz_aggr + 2;
    subsizes[0] = nx + 2;
    subsizes[1] = ny + 2;
    subsizes[2] = nz + 2;
    starts[0] = 0;
    starts[1] = 0;
    starts[2] = 0;
    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &ddtype_temp1);
    MPI_Type_size(MPI_DOUBLE, &r8size);
    lb = 0;
    extent = (MPI_Aint) r8size;
    MPI_Type_create_resized(ddtype_temp1, lb, extent, &ddtype_scatterv);
    MPI_Type_commit(&ddtype_scatterv);

    starts[0] = ix;
    starts[1] = iy;
    starts[2] = iz;

    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &ddtype_temp1);
    MPI_Type_size(MPI_DOUBLE, &r8size);
    lb = 0;
    extent = (MPI_Aint) r8size;
    MPI_Type_create_resized(ddtype_temp1, lb, extent, &ddtype_recv);
    MPI_Type_commit(&ddtype_recv);

    for (i = 0; i < nprocs; i++) {
        cnt_scatterv[i] = 1;
        disps_scatterv[i] = nz * (cart_coord[i][2])
                          + ny * (cart_coord[i][1]) * (nz_aggr + 2)       
                          + nx * (cart_coord[i][0]) * (nz_aggr + 2) * (ny_aggr + 2);                            // x offset
    }

    nx = sdm->nx;
    ny = sdm->ny;
    nz = sdm->nz;
    size_t size = (nx+2)*(ny+2)*(nz+2);
    double *rsd = (double*) calloc(size, sizeof(double));

    multigrid_residual(rsd, a_poisson, sol, rhs, sdm, sdm->is_aggregated);
    vv_dot_3d_matrix(&res0tol, rsd, rsd, nx, ny, nz, sdm->is_aggregated);

    for(cyc = 1; cyc <= n_vcycles; cyc++)
    {
        rbgs_iterator_poisson_matrix(sol, a_poisson, rhs, sdm, maxiteration, omega_sor, sdm->is_aggregated);
        multigrid_residual(rsd, a_poisson, sol, rhs, sdm, sdm->is_aggregated);
        multigrid_restriction(mg_sdm[1].b, rsd, &mg_sdm[1], sdm, 0);
        if(myrank==0) printf("[MG] Restriction from level 0 to level 1\n");

        for(l=1; l<=lv_aggregation-1; l++)
        {
            size_t level_size = (mg_sdm[l].nx + 2) * (mg_sdm[l].ny + 2) * (mg_sdm[l].nz + 2);
            memset(mg_sdm[l].x, 0, level_size * sizeof(double));
            rbgs_iterator_poisson_matrix(mg_sdm[l].x, &mg_a_poisson[l], mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
            multigrid_residual(mg_sdm[l].r, &mg_a_poisson[l], mg_sdm[l].x, mg_sdm[l].b, &mg_sdm[l], mg_sdm[l].is_aggregated);
            multigrid_restriction(mg_sdm[l+1].b, mg_sdm[l].r, &mg_sdm[l+1], &mg_sdm[l], l);
            if(myrank==0) printf("[MG] Restriction from level %d to level %d\n", l, l+1);
        }
    
l = lv_aggregation;

if (myrank == 0) {
    MPI_Gatherv(MPI_IN_PLACE, 0, ddtype_send, mg_sdm[lv_aggregation].b, cnt_gatherv, disps_gatherv, ddtype_gatherv, 0, MPI_COMM_WORLD);
} else {
    MPI_Gatherv(mg_sdm[lv_aggregation].b, 1, ddtype_send, mg_sdm[lv_aggregation].b, cnt_gatherv, disps_gatherv, ddtype_gatherv, 0, MPI_COMM_WORLD);
}

    if (myrank == 0) {
    // Restriction phase
    l = lv_aggregation;

    for (int l = lv_aggregation; l < n_levels; l++) {
        size_t level_size = (mg_sdm[l].nx + 2) * (mg_sdm[l].ny + 2) * (mg_sdm[l].nz + 2);
        memset(mg_sdm[l].x, 0, level_size * sizeof(double));
        rbgs_iterator_poisson_matrix(mg_sdm[l].x, &mg_a_poisson[l], mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
        multigrid_residual(mg_sdm[l].r, &mg_a_poisson[l], mg_sdm[l].x, mg_sdm[l].b, &mg_sdm[l], mg_sdm[l].is_aggregated);
        multigrid_restriction(mg_sdm[l+1].b, mg_sdm[l].r, &mg_sdm[l+1], &mg_sdm[l], l);
        printf("[MG] Restriction from level %d to level %d\n", l, l+1);
    }

    // Solve at coarsest level
    multigrid_solve_coarset_level(mg_sdm[n_levels].x, &mg_a_poisson[n_levels], mg_sdm[n_levels].b, &mg_sdm[n_levels], 1000, tolerance, omega_sor, mg_sdm[n_levels].is_aggregated);

    multigrid_residual(mg_sdm[n_levels].r, &mg_a_poisson[n_levels], mg_sdm[n_levels].x, mg_sdm[n_levels].b, &mg_sdm[n_levels], mg_sdm[n_levels].is_aggregated);

    printf("[MG] Solution in the coarsest level : x(1,1,1) = %.10e, residue = %.10e\n", mg_sdm[n_levels].x[IDX(1,1,1,mg_sdm[n_levels].nx,mg_sdm[n_levels].ny,mg_sdm[n_levels].nz)], mg_sdm[n_levels].r[IDX(1,1,1,mg_sdm[n_levels].nx,mg_sdm[n_levels].ny,mg_sdm[n_levels].nz)]);

#ifdef DEBUG_COARSEST
    multigrid_common_print_coarsest_level_solution(cyc, &mg_sdm[n_levels]);
#endif

    // Prolongation phase
    for (l = n_levels-1; l >= lv_aggregation; l--) {
        multigrid_prolongation_linear_on_nonuniform_grid(mg_sdm[l].r, mg_sdm[l+1].x, &mg_sdm[l], &mg_sdm[l+1], l);

        // x = x + r
        nx = mg_sdm[l].nx;
        ny = mg_sdm[l].ny;
        nz = mg_sdm[l].nz;
        for (i = 0; i < (nx+2)*(ny+2)*(nz+2); i++) {
            mg_sdm[l].x[i] += mg_sdm[l].r[i];
        }

        rbgs_iterator_poisson_matrix(mg_sdm[l].x, &mg_a_poisson[l], mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
    }
}


if (myrank == 0) {
    MPI_Scatterv(mg_sdm[lv_aggregation].x, cnt_scatterv, disps_scatterv, ddtype_scatterv, MPI_IN_PLACE, 0, ddtype_recv, 0, MPI_COMM_WORLD);
} else {
    MPI_Scatterv(mg_sdm[lv_aggregation].x, cnt_scatterv, disps_scatterv, ddtype_scatterv, mg_sdm[lv_aggregation].x, 1, ddtype_recv, 0, MPI_COMM_WORLD);
}


// Loop for prolongation from aggregation level down to level 1
for (int l = lv_aggregation-1; l >= 1; l--) {
    geometry_halocell_update_selectively(mg_sdm[l+1].x, &mg_sdm[l+1], mg_sdm[l+1].is_aggregated);

    multigrid_prolongation_linear_on_nonuniform_grid(mg_sdm[l].r, mg_sdm[l+1].x,
                                                     &mg_sdm[l], &mg_sdm[l+1], l);

    nx = mg_sdm[l].nx;
    ny = mg_sdm[l].ny;
    nz = mg_sdm[l].nz;
    int size = (nx+2)*(ny+2)*(nz+2);
    for (i = 0; i < size; i++) {
        mg_sdm[l].x[i] += mg_sdm[l].r[i];
    }

    rbgs_iterator_poisson_matrix(mg_sdm[l].x, &mg_a_poisson[l], mg_sdm[l].b,
                                 &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
}

// Prolongation to actual solution
geometry_halocell_update_selectively(mg_sdm[1].x, &mg_sdm[1], mg_sdm[1].is_aggregated);

multigrid_prolongation_linear_on_nonuniform_grid(rsd, mg_sdm[1].x, sdm, &mg_sdm[1], 0);

nx = sdm->nx;
ny = sdm->ny;
nz = sdm->nz;
int size = (nx+2)*(ny+2)*(nz+2);
for (i = 0; i < size; i++) {
    sol[i] += rsd[i];
}

rbgs_iterator_poisson_matrix(sol, a_poisson, rhs, sdm, maxiteration, omega_sor, sdm->is_aggregated);

multigrid_residual(rsd, a_poisson, sol, rhs, sdm, sdm->is_aggregated);

vv_dot_3d_matrix(&rsd_val, rsd, rsd, nx, ny, nz, sdm->is_aggregated);

if (myrank == 0) {
    printf("[MG] cycle = %d, Error: %e %e %e\n", cyc, sqrt(rsd_val), sqrt(res0tol), sqrt(rsd_val/res0tol));
}

if (sqrt(rsd_val/res0tol) < tolerance) break;

    }

    if(myrank==0) printf("[MG] Total %d V-cycles end\n", cyc);

    free(rsd);
    rsd = NULL;

}

void multigrid_adaptive_aggregation_vcycle_solver(double *sol, matrix_poisson *a_poisson, double *rhs, subdomain *sdm, int maxiteration, double tolerance, double omega_sor)
{
    int l, i, cyc;
    double rsd_val, res0tol;

    MPI_Datatype  ddtype_temp1;
    MPI_Datatype ddtype_send_x, ddtype_send_y, ddtype_send_z;
    MPI_Datatype ddtype_gatherv_x, ddtype_gatherv_y, ddtype_gatherv_z;

    MPI_Datatype *ddtype_send_max_ptr;
    MPI_Datatype *ddtype_send_med_ptr;
    MPI_Datatype *ddtype_send_min_ptr;

    MPI_Datatype *ddtype_gatherv_max_ptr;
    MPI_Datatype *ddtype_gatherv_med_ptr;
    MPI_Datatype *ddtype_gatherv_min_ptr;

    int sizes[3], subsizes[3], starts[3];
    int r8size;
    int nx_aggr, ny_aggr, nz_aggr;
    int nx, ny, nz;
    int ix, iy, iz;
    int level_case;
    char max, med, min;

    MPI_Aint extent, lb;

    cart_comm_1d *comm_max_ptr = NULL, *comm_med_ptr = NULL, *comm_min_ptr = NULL;


    int *lv_aggr_max_ptr = NULL, *lv_aggr_med_ptr = NULL, *lv_aggr_min_ptr = NULL;
    
    int *cnt_gatherv_x  = NULL, *disps_gatherv_x  = NULL;
    int *cnt_gatherv_y  = NULL, *disps_gatherv_y  = NULL;
    int *cnt_gatherv_z  = NULL, *disps_gatherv_z  = NULL;

    int *cnt_gatherv_max_ptr   = NULL, *disps_gatherv_max_ptr   = NULL;
    int *cnt_gatherv_med_ptr   = NULL, *disps_gatherv_med_ptr   = NULL;
    int *cnt_gatherv_min_ptr   = NULL, *disps_gatherv_min_ptr   = NULL;

    int (*cart_coord)[3];
    cart_coord = malloc(nprocs * sizeof(*cart_coord));
    int ierr;
    for (i = 0; i < nprocs; i++) {
        ierr = MPI_Cart_coords(mpi_world_cart, i, 3, cart_coord[i]);
        if (ierr != MPI_SUCCESS) {
            fprintf(stderr, "MPI_Cart_coords error\n");
            MPI_Abort(MPI_COMM_WORLD, ierr);
        }
    }
    cnt_gatherv_x = malloc((comm_1d_x.nprocs ) * sizeof(int));
    disps_gatherv_x = malloc((comm_1d_x.nprocs ) * sizeof(int));
    cnt_gatherv_y = malloc((comm_1d_y.nprocs ) * sizeof(int));
    disps_gatherv_y = malloc((comm_1d_y.nprocs ) * sizeof(int));
    cnt_gatherv_z = malloc((comm_1d_z.nprocs ) * sizeof(int));
    disps_gatherv_z = malloc((comm_1d_z.nprocs ) * sizeof(int));



    if (lv_aggregation_x != 0)
    {
        nx_aggr = mg_sdm[lv_aggregation_x].nx;
        ny_aggr = mg_sdm[lv_aggregation_x].ny;
        nz_aggr = mg_sdm[lv_aggregation_x].nz;
    
        nx = nx_aggr / comm_1d_x.nprocs;
        ny = ny_aggr;
        nz = nz_aggr;

        ix = nx * comm_1d_x.myrank;
        iy = 0;
        iz = 0;

        // For MPI_Gatherv
        sizes[0] = nx_aggr + 2;
        sizes[1] = ny_aggr + 2;
        sizes[2] = nz_aggr + 2;
        subsizes[0] = nx;
        subsizes[1] = ny;
        subsizes[2] = nz;
        starts[0] = ix+1;
        starts[1] = iy+1;
        starts[2] = iz+1;
        MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &ddtype_temp1);
        MPI_Type_size(MPI_DOUBLE, &r8size);
        lb = 0;
        extent = (MPI_Aint) r8size;
        MPI_Type_create_resized(ddtype_temp1, lb, extent, &ddtype_send_x);
        MPI_Type_commit(&ddtype_send_x);
    
        starts[0] = 1;
        starts[1] = 1;
        starts[2] = 1;
        MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &ddtype_temp1);
        MPI_Type_size(MPI_DOUBLE, &r8size);
        lb = 0;
        extent = (MPI_Aint) r8size;
        MPI_Type_create_resized(ddtype_temp1, lb, extent, &ddtype_gatherv_x);
        MPI_Type_commit(&ddtype_gatherv_x);

        for (i = 0; i < comm_1d_x.nprocs; i++) {
            cnt_gatherv_x[i] = 1;
            disps_gatherv_x[i] = nx * i * (ny+2) * (nz+2);    
        }
    }

    if (lv_aggregation_y != 0)
    {
        nx_aggr = mg_sdm[lv_aggregation_y].nx;
        ny_aggr = mg_sdm[lv_aggregation_y].ny;
        nz_aggr = mg_sdm[lv_aggregation_y].nz;
    
        nx = nx_aggr;
        ny = ny_aggr / comm_1d_y.nprocs;
        nz = nz_aggr;

        ix = 0;
        iy = ny * comm_1d_y.myrank;
        iz = 0;

        // For MPI_Gatherv
        sizes[0] = nx_aggr + 2;
        sizes[1] = ny_aggr + 2;
        sizes[2] = nz_aggr + 2;
        subsizes[0] = nx;
        subsizes[1] = ny;
        subsizes[2] = nz;
        starts[0] = ix+1;
        starts[1] = iy+1;
        starts[2] = iz+1;
        MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &ddtype_temp1);
        MPI_Type_size(MPI_DOUBLE, &r8size);
        lb = 0;
        extent = (MPI_Aint) r8size;
        MPI_Type_create_resized(ddtype_temp1, lb, extent, &ddtype_send_y);
        MPI_Type_commit(&ddtype_send_y);
    
        starts[0] = 1;
        starts[1] = 1;
        starts[2] = 1;
        MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &ddtype_temp1);
        MPI_Type_size(MPI_DOUBLE, &r8size);
        lb = 0;
        extent = (MPI_Aint) r8size;
        MPI_Type_create_resized(ddtype_temp1, lb, extent, &ddtype_gatherv_y);
        MPI_Type_commit(&ddtype_gatherv_y);

        for (i = 0; i < comm_1d_y.nprocs; i++) {
            cnt_gatherv_y[i] = 1;
            disps_gatherv_y[i] = ny * i * (nz+2);    
        }
    }

    if (lv_aggregation_z != 0)
    {
        nx_aggr = mg_sdm[lv_aggregation_z].nx;
        ny_aggr = mg_sdm[lv_aggregation_z].ny;
        nz_aggr = mg_sdm[lv_aggregation_z].nz;
    
        nx = nx_aggr;
        ny = ny_aggr;
        nz = nz_aggr / comm_1d_z.nprocs;

        ix = 0;
        iy = 0;
        iz = nz * comm_1d_z.myrank;

        // For MPI_Gatherv
        sizes[0] = nx_aggr + 2;
        sizes[1] = ny_aggr + 2;
        sizes[2] = nz_aggr + 2;
        subsizes[0] = nx;
        subsizes[1] = ny;
        subsizes[2] = nz;
        starts[0] = ix+1;
        starts[1] = iy+1;
        starts[2] = iz+1;
        MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &ddtype_temp1);
        MPI_Type_size(MPI_DOUBLE, &r8size);
        lb = 0;
        extent = (MPI_Aint) r8size;
        MPI_Type_create_resized(ddtype_temp1, lb, extent, &ddtype_send_z);
        MPI_Type_commit(&ddtype_send_z);
    
        starts[0] = 1;
        starts[1] = 1;
        starts[2] = 1;
        MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &ddtype_temp1);
        MPI_Type_size(MPI_DOUBLE, &r8size);
        lb = 0;
        extent = (MPI_Aint) r8size;
        MPI_Type_create_resized(ddtype_temp1, lb, extent, &ddtype_gatherv_z);
        MPI_Type_commit(&ddtype_gatherv_z);

        for (i = 0; i < comm_1d_z.nprocs; i++) {
            cnt_gatherv_z[i] = 1;
            disps_gatherv_z[i] = nz * i;    
        }
    }

    if (lv_aggregation_x == 0)
    {
        if (lv_aggregation_y == 0)
        {
            level_case = 1;
            lv_aggr_max_ptr = &lv_aggregation_z;
            comm_max_ptr = &comm_1d_z;
            ddtype_send_max_ptr = &ddtype_send_z;
            ddtype_gatherv_max_ptr = &ddtype_gatherv_z;
            cnt_gatherv_max_ptr = cnt_gatherv_z;
            disps_gatherv_max_ptr = disps_gatherv_z;
            max = 'z';
        }
        else if ( lv_aggregation_z == 0 )
        {
            level_case = 1;
            lv_aggr_max_ptr = &lv_aggregation_y;
            comm_max_ptr = &comm_1d_y;
            ddtype_send_max_ptr = &ddtype_send_y;
            ddtype_gatherv_max_ptr = &ddtype_gatherv_y;
            cnt_gatherv_max_ptr = cnt_gatherv_y;
            disps_gatherv_max_ptr = disps_gatherv_y;
            max = 'y';
        }
        else
        {
            level_case = 2;
            if (lv_aggregation_y > lv_aggregation_z)
            {
                lv_aggr_max_ptr = &lv_aggregation_y;
                comm_max_ptr = &comm_1d_y;
                ddtype_send_max_ptr = &ddtype_send_y;
                ddtype_gatherv_max_ptr = &ddtype_gatherv_y;
                cnt_gatherv_max_ptr = cnt_gatherv_y;
                disps_gatherv_max_ptr = disps_gatherv_y;
                max = 'y';
                lv_aggr_med_ptr = &lv_aggregation_z;
                comm_med_ptr = &comm_1d_z;
                ddtype_send_med_ptr = &ddtype_send_z;
                ddtype_gatherv_med_ptr = &ddtype_gatherv_z;
                cnt_gatherv_med_ptr = cnt_gatherv_z;
                disps_gatherv_med_ptr = disps_gatherv_z;
                med = 'z';
            }
            else
            {
                lv_aggr_max_ptr = &lv_aggregation_z;
                comm_max_ptr = &comm_1d_z;
                ddtype_send_max_ptr = &ddtype_send_z;
                ddtype_gatherv_max_ptr = &ddtype_gatherv_z;
                cnt_gatherv_max_ptr = cnt_gatherv_z;
                disps_gatherv_max_ptr = disps_gatherv_z;
                max = 'z';
                lv_aggr_med_ptr = &lv_aggregation_y;
                comm_med_ptr = &comm_1d_y;
                ddtype_send_med_ptr = &ddtype_send_y;
                ddtype_gatherv_med_ptr = &ddtype_gatherv_y;
                cnt_gatherv_med_ptr = cnt_gatherv_y;
                disps_gatherv_med_ptr = disps_gatherv_y;
                med = 'y';
            }
        }
    }

    else if ( lv_aggregation_y == 0 )
    {
        if (lv_aggregation_z == 0)
        {
            level_case = 1;
            lv_aggr_max_ptr = &lv_aggregation_x;
            comm_max_ptr = &comm_1d_x;
            ddtype_send_max_ptr = &ddtype_send_x;
            ddtype_gatherv_max_ptr = &ddtype_gatherv_x;
            cnt_gatherv_max_ptr = cnt_gatherv_x;
            disps_gatherv_max_ptr = disps_gatherv_x;
            max = 'x';
        }
        else
        {
            level_case = 2;
            if (lv_aggregation_x > lv_aggregation_z)
            {
                lv_aggr_max_ptr = &lv_aggregation_x;
                comm_max_ptr = &comm_1d_x;
                ddtype_send_max_ptr = &ddtype_send_x;
                ddtype_gatherv_max_ptr = &ddtype_gatherv_x;
                cnt_gatherv_max_ptr = cnt_gatherv_x;
                disps_gatherv_max_ptr = disps_gatherv_x;
                max = 'x';
                lv_aggr_med_ptr = &lv_aggregation_z;
                comm_med_ptr = &comm_1d_z;
                ddtype_send_med_ptr = &ddtype_send_z;
                ddtype_gatherv_med_ptr = &ddtype_gatherv_z;
                cnt_gatherv_med_ptr = cnt_gatherv_z;
                disps_gatherv_med_ptr = disps_gatherv_z;
                med = 'z';
            }
            else
            {
                lv_aggr_max_ptr = &lv_aggregation_z;
                comm_max_ptr = &comm_1d_z;
                ddtype_send_max_ptr = &ddtype_send_z;
                ddtype_gatherv_max_ptr = &ddtype_gatherv_z;
                cnt_gatherv_max_ptr = cnt_gatherv_z;
                disps_gatherv_max_ptr = disps_gatherv_z;
                max = 'z';
                lv_aggr_med_ptr = &lv_aggregation_x;
                comm_med_ptr = &comm_1d_x;
                ddtype_send_med_ptr = &ddtype_send_x;
                ddtype_gatherv_med_ptr = &ddtype_gatherv_x;
                cnt_gatherv_med_ptr = cnt_gatherv_x;
                disps_gatherv_med_ptr = disps_gatherv_x;
                med = 'x';
            }
        }
    }

    else if ( lv_aggregation_z == 0 )
    {
        level_case = 2;
        if (lv_aggregation_x > lv_aggregation_y)
        {
            lv_aggr_max_ptr = &lv_aggregation_x;
            comm_max_ptr = &comm_1d_x;
            ddtype_send_max_ptr = &ddtype_send_x;
            ddtype_gatherv_max_ptr = &ddtype_gatherv_x;
            cnt_gatherv_max_ptr = cnt_gatherv_x;
            disps_gatherv_max_ptr = disps_gatherv_x;
            max = 'x';
            lv_aggr_med_ptr = &lv_aggregation_y;
            comm_med_ptr = &comm_1d_y;
            ddtype_send_med_ptr = &ddtype_send_y;
            ddtype_gatherv_med_ptr = &ddtype_gatherv_y;
            cnt_gatherv_med_ptr = cnt_gatherv_y;
            disps_gatherv_med_ptr = disps_gatherv_y;
            med = 'y';
        }
        else
        {
            lv_aggr_max_ptr = &lv_aggregation_y;
            comm_max_ptr = &comm_1d_y;
            ddtype_send_max_ptr = &ddtype_send_y;
            ddtype_gatherv_max_ptr = &ddtype_gatherv_y;
            cnt_gatherv_max_ptr = cnt_gatherv_y;
            disps_gatherv_max_ptr = disps_gatherv_y;
            max = 'y';
            lv_aggr_med_ptr = &lv_aggregation_x;
            comm_med_ptr = &comm_1d_x;
            ddtype_send_med_ptr = &ddtype_send_x;
            ddtype_gatherv_med_ptr = &ddtype_gatherv_x;
            cnt_gatherv_med_ptr = cnt_gatherv_x;
            disps_gatherv_med_ptr = disps_gatherv_x;
            med = 'x';
        }
    }

    else
    {
        level_case = 3;

        if (lv_aggregation_x > lv_aggregation_y)
        {
            if (lv_aggregation_x > lv_aggregation_z)
            {
                if (lv_aggregation_y > lv_aggregation_z)
                {
                    /* x > y > z */
                    lv_aggr_max_ptr = &lv_aggregation_x;
                    lv_aggr_med_ptr = &lv_aggregation_y;
                    lv_aggr_min_ptr = &lv_aggregation_z;
                    comm_max_ptr = &comm_1d_x;
                    comm_med_ptr = &comm_1d_y;
                    comm_min_ptr = &comm_1d_z;
                    ddtype_send_max_ptr    = &ddtype_send_x;
                    ddtype_gatherv_max_ptr = &ddtype_gatherv_x;
                    cnt_gatherv_max_ptr    = cnt_gatherv_x;
                    disps_gatherv_max_ptr  = disps_gatherv_x;
                    ddtype_send_med_ptr    = &ddtype_send_y;
                    ddtype_gatherv_med_ptr = &ddtype_gatherv_y;
                    cnt_gatherv_med_ptr    = cnt_gatherv_y;
                    disps_gatherv_med_ptr  = disps_gatherv_y;
                    ddtype_send_min_ptr    = &ddtype_send_z;
                    ddtype_gatherv_min_ptr = &ddtype_gatherv_z;
                    cnt_gatherv_min_ptr    = cnt_gatherv_z;
                    disps_gatherv_min_ptr  = disps_gatherv_z;
                    max = 'x'; med = 'y'; min = 'z';
                }
                else
                {
                    /* x > z >= y */
                    lv_aggr_max_ptr = &lv_aggregation_x;
                    lv_aggr_med_ptr = &lv_aggregation_z;
                    lv_aggr_min_ptr = &lv_aggregation_y;
                    comm_max_ptr = &comm_1d_x;
                    comm_med_ptr = &comm_1d_z;
                    comm_min_ptr = &comm_1d_y;
                    ddtype_send_max_ptr    = &ddtype_send_x;
                    ddtype_gatherv_max_ptr = &ddtype_gatherv_x;
                    cnt_gatherv_max_ptr    = cnt_gatherv_x;
                    disps_gatherv_max_ptr  = disps_gatherv_x;
                    ddtype_send_med_ptr    = &ddtype_send_z;
                    ddtype_gatherv_med_ptr = &ddtype_gatherv_z;
                    cnt_gatherv_med_ptr    = cnt_gatherv_z;
                    disps_gatherv_med_ptr  = disps_gatherv_z;
                    ddtype_send_min_ptr    = &ddtype_send_y;
                    ddtype_gatherv_min_ptr = &ddtype_gatherv_y;
                    cnt_gatherv_min_ptr    = cnt_gatherv_y;
                    disps_gatherv_min_ptr  = disps_gatherv_y;
                    max = 'x'; med = 'z'; min = 'y';
                }
            }
            else
            {
                /* z >= x > y */
                lv_aggr_max_ptr = &lv_aggregation_z;
                lv_aggr_med_ptr = &lv_aggregation_x;
                lv_aggr_min_ptr = &lv_aggregation_y;
                comm_max_ptr = &comm_1d_z;
                comm_med_ptr = &comm_1d_x;
                comm_min_ptr = &comm_1d_y;
                ddtype_send_max_ptr    = &ddtype_send_z;
                ddtype_gatherv_max_ptr = &ddtype_gatherv_z;
                cnt_gatherv_max_ptr    = cnt_gatherv_z;
                disps_gatherv_max_ptr  = disps_gatherv_z;
                ddtype_send_med_ptr    = &ddtype_send_x;
                ddtype_gatherv_med_ptr = &ddtype_gatherv_x;
                cnt_gatherv_med_ptr    = cnt_gatherv_x;
                disps_gatherv_med_ptr  = disps_gatherv_x;
                ddtype_send_min_ptr    = &ddtype_send_y;
                ddtype_gatherv_min_ptr = &ddtype_gatherv_y;
                cnt_gatherv_min_ptr    = cnt_gatherv_y;
                disps_gatherv_min_ptr  = disps_gatherv_y;
                max = 'z'; med = 'x'; min = 'y';
            }
        }
        else
        {
            if (lv_aggregation_y > lv_aggregation_z)
            {
                if (lv_aggregation_x > lv_aggregation_z)
                {
                    /* y >= x > z */
                    lv_aggr_max_ptr = &lv_aggregation_y;
                    lv_aggr_med_ptr = &lv_aggregation_x;
                    lv_aggr_min_ptr = &lv_aggregation_z;
                    comm_max_ptr = &comm_1d_y;
                    comm_med_ptr = &comm_1d_x;
                    comm_min_ptr = &comm_1d_z;
                    ddtype_send_max_ptr    = &ddtype_send_y;
                    ddtype_gatherv_max_ptr = &ddtype_gatherv_y;
                    cnt_gatherv_max_ptr    = cnt_gatherv_y;
                    disps_gatherv_max_ptr  = disps_gatherv_y;
                    ddtype_send_med_ptr    = &ddtype_send_x;
                    ddtype_gatherv_med_ptr = &ddtype_gatherv_x;
                    cnt_gatherv_med_ptr    = cnt_gatherv_x;
                    disps_gatherv_med_ptr  = disps_gatherv_x;
                    ddtype_send_min_ptr    = &ddtype_send_z;
                    ddtype_gatherv_min_ptr = &ddtype_gatherv_z;
                    cnt_gatherv_min_ptr    = cnt_gatherv_z;
                    disps_gatherv_min_ptr  = disps_gatherv_z;
                    max = 'y'; med = 'x'; min = 'z';
                }
                else
                {
                    /* y >= z >= x */
                    lv_aggr_max_ptr = &lv_aggregation_y;
                    lv_aggr_med_ptr = &lv_aggregation_z;
                    lv_aggr_min_ptr = &lv_aggregation_x;
                    comm_max_ptr = &comm_1d_y;
                    comm_med_ptr = &comm_1d_z;
                    comm_min_ptr = &comm_1d_x;
                    ddtype_send_max_ptr    = &ddtype_send_y;
                    ddtype_gatherv_max_ptr = &ddtype_gatherv_y;
                    cnt_gatherv_max_ptr    = cnt_gatherv_y;
                    disps_gatherv_max_ptr  = disps_gatherv_y;
                    ddtype_send_med_ptr    = &ddtype_send_z;
                    ddtype_gatherv_med_ptr = &ddtype_gatherv_z;
                    cnt_gatherv_med_ptr    = cnt_gatherv_z;
                    disps_gatherv_med_ptr  = disps_gatherv_z;
                    ddtype_send_min_ptr    = &ddtype_send_x;
                    ddtype_gatherv_min_ptr = &ddtype_gatherv_x;
                    cnt_gatherv_min_ptr    = cnt_gatherv_x;
                    disps_gatherv_min_ptr  = disps_gatherv_x;
                    max = 'y'; med = 'z'; min = 'x';
                }
            }
            else
            {
                /* z >= y >= x */
                lv_aggr_max_ptr = &lv_aggregation_z;
                lv_aggr_med_ptr = &lv_aggregation_y;
                lv_aggr_min_ptr = &lv_aggregation_x;
                comm_max_ptr = &comm_1d_z;
                comm_med_ptr = &comm_1d_y;
                comm_min_ptr = &comm_1d_x;
                ddtype_send_max_ptr    = &ddtype_send_z;
                ddtype_gatherv_max_ptr = &ddtype_gatherv_z;
                cnt_gatherv_max_ptr    = cnt_gatherv_z;
                disps_gatherv_max_ptr  = disps_gatherv_z;
                ddtype_send_med_ptr    = &ddtype_send_y;
                ddtype_gatherv_med_ptr = &ddtype_gatherv_y;
                cnt_gatherv_med_ptr    = cnt_gatherv_y;
                disps_gatherv_med_ptr  = disps_gatherv_y;
                ddtype_send_min_ptr    = &ddtype_send_x;
                ddtype_gatherv_min_ptr = &ddtype_gatherv_x;
                cnt_gatherv_min_ptr    = cnt_gatherv_x;
                disps_gatherv_min_ptr  = disps_gatherv_x;
                max = 'z'; med = 'y'; min = 'x';
            }
        }
        if(myrank==0) printf("myrank=%d, lv_aggr_max_ptr= %d, lv_aggr_med_ptr= %d, lv_aggr_min_ptr= %d, max=%c, med=%c, min=%c\n",myrank, *lv_aggr_max_ptr, *lv_aggr_med_ptr, *lv_aggr_min_ptr, max, med, min );
    }
    if(myrank==0) printf("level_case=%d\n", level_case);

    nx = sdm->nx;
    ny = sdm->ny;
    nz = sdm->nz;
    size_t size = (nx+2)*(ny+2)*(nz+2);
    double *rsd = (double*) calloc(size, sizeof(double));

    switch (level_case) {
        case 1:
            // 1 aggregation level

            multigrid_residual(rsd, a_poisson, sol, rhs, sdm, sdm->is_aggregated);
            vv_dot_3d_matrix(&res0tol, rsd, rsd, nx, ny, nz, sdm->is_aggregated);

            for(cyc = 1; cyc <= n_vcycles; cyc++)
            {
                rbgs_iterator_poisson_matrix(sol, a_poisson, rhs, sdm, maxiteration, omega_sor, sdm->is_aggregated);
                multigrid_residual(rsd, a_poisson, sol, rhs, sdm, sdm->is_aggregated);
                multigrid_restriction(mg_sdm[1].b, rsd, &mg_sdm[1], sdm, 0);
                if(myrank==0) printf("[MG] Restriction from level 0 to level 1\n");

                for(l=1; l<=*lv_aggr_max_ptr-1; l++)
                {
                    size_t level_size = (mg_sdm[l].nx + 2) * (mg_sdm[l].ny + 2) * (mg_sdm[l].nz + 2);
                    memset(mg_sdm[l].x, 0, level_size * sizeof(double));
                    rbgs_iterator_poisson_matrix(mg_sdm[l].x, &mg_a_poisson[l], mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                    multigrid_residual(mg_sdm[l].r, &mg_a_poisson[l], mg_sdm[l].x, mg_sdm[l].b, &mg_sdm[l], mg_sdm[l].is_aggregated);
                    multigrid_restriction(mg_sdm[l+1].b, mg_sdm[l].r, &mg_sdm[l+1], &mg_sdm[l], l);
                    if(myrank==0) printf("[MG] Restriction from level %d to level %d\n", l, l+1);
                }
    
l = *lv_aggr_max_ptr;

if (myrank == 0) {
    MPI_Gatherv(MPI_IN_PLACE, 0, *ddtype_send_max_ptr, mg_sdm[l].b, cnt_gatherv_max_ptr, disps_gatherv_max_ptr, *ddtype_gatherv_max_ptr, 0, comm_max_ptr->mpi_comm);
} else {
    MPI_Gatherv(mg_sdm[l].b, 1, *ddtype_send_max_ptr, mg_sdm[l].b, cnt_gatherv_max_ptr, disps_gatherv_max_ptr, *ddtype_gatherv_max_ptr, 0, comm_max_ptr->mpi_comm);
}

                if ((*comm_max_ptr).myrank == 0) {
                    // Restriction phase
                    for (int l = *lv_aggr_max_ptr; l < n_levels; l++) {
                        size_t level_size = (mg_sdm[l].nx + 2) * (mg_sdm[l].ny + 2) * (mg_sdm[l].nz + 2);
                        memset(mg_sdm[l].x, 0, level_size * sizeof(double));
                        rbgs_iterator_poisson_matrix(mg_sdm[l].x, &mg_a_poisson[l], mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                        multigrid_residual(mg_sdm[l].r, &mg_a_poisson[l], mg_sdm[l].x, mg_sdm[l].b, &mg_sdm[l], mg_sdm[l].is_aggregated);
                        multigrid_restriction(mg_sdm[l+1].b, mg_sdm[l].r, &mg_sdm[l+1], &mg_sdm[l], l);
                        printf("[MG] Restriction from level %d to level %d\n", l, l+1);
                    }

                    // Solve at coarsest level
                    multigrid_solve_coarset_level(mg_sdm[n_levels].x, &mg_a_poisson[n_levels], mg_sdm[n_levels].b, &mg_sdm[n_levels], 1000, tolerance, omega_sor, mg_sdm[n_levels].is_aggregated);

                    multigrid_residual(mg_sdm[n_levels].r, &mg_a_poisson[n_levels], mg_sdm[n_levels].x, mg_sdm[n_levels].b, &mg_sdm[n_levels], mg_sdm[n_levels].is_aggregated);

                    printf("[MG] Solution in the coarsest level : x(1,1,1) = %.10e, residue = %.10e\n", mg_sdm[n_levels].x[IDX(1,1,1,mg_sdm[n_levels].nx,mg_sdm[n_levels].ny,mg_sdm[n_levels].nz)], mg_sdm[n_levels].r[IDX(1,1,1,mg_sdm[n_levels].nx,mg_sdm[n_levels].ny,mg_sdm[n_levels].nz)]);

#ifdef DEBUG_COARSEST
    multigrid_common_print_coarsest_level_solution(cyc, &mg_sdm[n_levels]);
#endif

                    // Prolongation phase
                    for (l = n_levels-1; l >= *lv_aggr_max_ptr; l--) {
                        multigrid_prolongation_linear_on_nonuniform_grid(mg_sdm[l].r, mg_sdm[l+1].x, &mg_sdm[l], &mg_sdm[l+1], l);

                        // x = x + r
                        nx = mg_sdm[l].nx;
                        ny = mg_sdm[l].ny;
                        nz = mg_sdm[l].nz;
                        for (i = 0; i < (nx+2)*(ny+2)*(nz+2); i++) {
                            mg_sdm[l].x[i] += mg_sdm[l].r[i];
                        }

                    rbgs_iterator_poisson_matrix(mg_sdm[l].x, &mg_a_poisson[l], mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                    }
                }

l = *lv_aggr_max_ptr;
// #ifdef MPI_IN_PLACE
if (myrank == 0) {
    MPI_Scatterv(mg_sdm[l].x, cnt_gatherv_max_ptr, disps_gatherv_max_ptr, *ddtype_gatherv_max_ptr, MPI_IN_PLACE, 0, *ddtype_send_max_ptr, 0, comm_max_ptr->mpi_comm);
} else {
    MPI_Scatterv(mg_sdm[l].x, cnt_gatherv_max_ptr, disps_gatherv_max_ptr, *ddtype_gatherv_max_ptr, mg_sdm[l].x, 1, *ddtype_send_max_ptr, 0, comm_max_ptr->mpi_comm);
}


                // Loop for prolongation from aggregation level down to level 1
                for (int l = *lv_aggr_max_ptr-1; l >= 1; l--) {
                    geometry_halocell_update_selectively(mg_sdm[l+1].x, &mg_sdm[l+1], mg_sdm[l+1].is_aggregated);

                    multigrid_prolongation_linear_on_nonuniform_grid(mg_sdm[l].r, mg_sdm[l+1].x,
                                                     &mg_sdm[l], &mg_sdm[l+1], l);
    

                    nx = mg_sdm[l].nx;
                    ny = mg_sdm[l].ny;
                    nz = mg_sdm[l].nz;
                    int size = (nx+2)*(ny+2)*(nz+2);
                    for (i = 0; i < size; i++) {
                        mg_sdm[l].x[i] += mg_sdm[l].r[i];
                    }

                    rbgs_iterator_poisson_matrix(mg_sdm[l].x, &mg_a_poisson[l], mg_sdm[l].b,
                                 &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                }

                // Prolongation to actual solution
                geometry_halocell_update_selectively(mg_sdm[1].x, &mg_sdm[1], mg_sdm[1].is_aggregated);

                multigrid_prolongation_linear_on_nonuniform_grid(rsd, mg_sdm[1].x, sdm, &mg_sdm[1], 0);

                nx = sdm->nx;
                ny = sdm->ny;
                nz = sdm->nz;
                int size = (nx+2)*(ny+2)*(nz+2);
                for (i = 0; i < size; i++) {
                    sol[i] += rsd[i];
                }

                rbgs_iterator_poisson_matrix(sol, a_poisson, rhs, sdm, maxiteration, omega_sor, sdm->is_aggregated);

                multigrid_residual(rsd, a_poisson, sol, rhs, sdm, sdm->is_aggregated);

                vv_dot_3d_matrix(&rsd_val, rsd, rsd, nx, ny, nz, sdm->is_aggregated);

                if (myrank == 0) {
                    printf("[MG] cycle = %d, Error: %e %e %e\n", cyc, sqrt(rsd_val), sqrt(res0tol), sqrt(rsd_val/res0tol));
                }

                if (sqrt(rsd_val/res0tol) < tolerance) break;


            }

            if(myrank==0) printf("[MG] Total %d V-cycles end\n", cyc);

            free(rsd);
            rsd = NULL;


            break;

        case 2:
            // 2 aggregation level

            multigrid_residual(rsd, a_poisson, sol, rhs, sdm, sdm->is_aggregated);
            vv_dot_3d_matrix(&res0tol, rsd, rsd, nx, ny, nz, sdm->is_aggregated);

            for(cyc = 1; cyc <= n_vcycles; cyc++)
            {
                rbgs_iterator_poisson_matrix(sol, a_poisson, rhs, sdm, maxiteration, omega_sor, sdm->is_aggregated);
                multigrid_residual(rsd, a_poisson, sol, rhs, sdm, sdm->is_aggregated);
                multigrid_restriction(mg_sdm[1].b, rsd, &mg_sdm[1], sdm, 0);
                if(myrank==0) printf("[MG] Restriction from level 0 to level 1\n");

                for(l=1; l<=*lv_aggr_med_ptr-1; l++)
                {
                    size_t level_size = (mg_sdm[l].nx + 2) * (mg_sdm[l].ny + 2) * (mg_sdm[l].nz + 2);
                    memset(mg_sdm[l].x, 0, level_size * sizeof(double));
                    rbgs_iterator_poisson_matrix(mg_sdm[l].x, &mg_a_poisson[l], mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                    multigrid_residual(mg_sdm[l].r, &mg_a_poisson[l], mg_sdm[l].x, mg_sdm[l].b, &mg_sdm[l], mg_sdm[l].is_aggregated);
                    multigrid_restriction(mg_sdm[l+1].b, mg_sdm[l].r, &mg_sdm[l+1], &mg_sdm[l], l);
                    if(myrank==0) printf("[MG] Restriction from level %d to level %d\n", l, l+1);
                }
    
l = *lv_aggr_med_ptr;

if (myrank == 0) {
    MPI_Gatherv(MPI_IN_PLACE, 0, *ddtype_send_med_ptr, mg_sdm[l].b, cnt_gatherv_med_ptr, disps_gatherv_med_ptr, *ddtype_gatherv_med_ptr, 0, comm_med_ptr->mpi_comm);
} else {
    MPI_Gatherv(mg_sdm[l].b, 1, *ddtype_send_med_ptr, mg_sdm[l].b, cnt_gatherv_med_ptr, disps_gatherv_med_ptr, *ddtype_gatherv_med_ptr, 0, comm_med_ptr->mpi_comm);
}

                if ((*comm_med_ptr).myrank == 0) {
                    // Restriction phase
                    for (int l = *lv_aggr_med_ptr; l < *lv_aggr_max_ptr; l++) {
                        size_t level_size = (mg_sdm[l].nx + 2) * (mg_sdm[l].ny + 2) * (mg_sdm[l].nz + 2);
                        memset(mg_sdm[l].x, 0, level_size * sizeof(double));
                        rbgs_iterator_poisson_matrix(mg_sdm[l].x, &mg_a_poisson[l], mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                        multigrid_residual(mg_sdm[l].r, &mg_a_poisson[l], mg_sdm[l].x, mg_sdm[l].b, &mg_sdm[l], mg_sdm[l].is_aggregated);
                        multigrid_restriction(mg_sdm[l+1].b, mg_sdm[l].r, &mg_sdm[l+1], &mg_sdm[l], l);
                        printf("[MG] Restriction from level %d to level %d\n", l, l+1);
                    }

l = *lv_aggr_max_ptr;

if (myrank == 0) {
    MPI_Gatherv(MPI_IN_PLACE, 0, *ddtype_send_max_ptr, mg_sdm[l].b, cnt_gatherv_max_ptr, disps_gatherv_max_ptr, *ddtype_gatherv_max_ptr, 0, comm_max_ptr->mpi_comm);
} else {
    MPI_Gatherv(mg_sdm[l].b, 1, *ddtype_send_max_ptr, mg_sdm[l].b, cnt_gatherv_max_ptr, disps_gatherv_max_ptr, *ddtype_gatherv_max_ptr, 0, comm_max_ptr->mpi_comm);
}

                if ((*comm_max_ptr).myrank == 0) {
                    // Restriction phase
                    for (int l = *lv_aggr_max_ptr; l < n_levels; l++) {
                        size_t level_size = (mg_sdm[l].nx + 2) * (mg_sdm[l].ny + 2) * (mg_sdm[l].nz + 2);
                        memset(mg_sdm[l].x, 0, level_size * sizeof(double));
                        rbgs_iterator_poisson_matrix(mg_sdm[l].x, &mg_a_poisson[l], mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                        multigrid_residual(mg_sdm[l].r, &mg_a_poisson[l], mg_sdm[l].x, mg_sdm[l].b, &mg_sdm[l], mg_sdm[l].is_aggregated);
                        multigrid_restriction(mg_sdm[l+1].b, mg_sdm[l].r, &mg_sdm[l+1], &mg_sdm[l], l);
                        printf("[MG] Restriction from level %d to level %d\n", l, l+1);
                    }

                    // Solve at coarsest level
                    multigrid_solve_coarset_level(mg_sdm[n_levels].x, &mg_a_poisson[n_levels], mg_sdm[n_levels].b, &mg_sdm[n_levels], 1000, tolerance, omega_sor, mg_sdm[n_levels].is_aggregated);

                    multigrid_residual(mg_sdm[n_levels].r, &mg_a_poisson[n_levels], mg_sdm[n_levels].x, mg_sdm[n_levels].b, &mg_sdm[n_levels], mg_sdm[n_levels].is_aggregated);

                    printf("[MG] Solution in the coarsest level : x(1,1,1) = %.10e, residue = %.10e\n", mg_sdm[n_levels].x[IDX(1,1,1,mg_sdm[n_levels].nx,mg_sdm[n_levels].ny,mg_sdm[n_levels].nz)], mg_sdm[n_levels].r[IDX(1,1,1,mg_sdm[n_levels].nx,mg_sdm[n_levels].ny,mg_sdm[n_levels].nz)]);

#ifdef DEBUG_COARSEST
    multigrid_common_print_coarsest_level_solution(cyc, &mg_sdm[n_levels]);
#endif

                    // Prolongation phase
                    for (l = n_levels-1; l >= *lv_aggr_max_ptr; l--) {
                        multigrid_prolongation_linear_on_nonuniform_grid(mg_sdm[l].r, mg_sdm[l+1].x, &mg_sdm[l], &mg_sdm[l+1], l);

                        // x = x + r
                        nx = mg_sdm[l].nx;
                        ny = mg_sdm[l].ny;
                        nz = mg_sdm[l].nz;
                        for (i = 0; i < (nx+2)*(ny+2)*(nz+2); i++) {
                            mg_sdm[l].x[i] += mg_sdm[l].r[i];
                        }

                    rbgs_iterator_poisson_matrix(mg_sdm[l].x, &mg_a_poisson[l], mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                    }
                }

l = *lv_aggr_max_ptr;
// #ifdef MPI_IN_PLACE
if (myrank == 0) {
    MPI_Scatterv(mg_sdm[l].x, cnt_gatherv_max_ptr, disps_gatherv_max_ptr, *ddtype_gatherv_max_ptr, MPI_IN_PLACE, 0, *ddtype_send_max_ptr, 0, comm_max_ptr->mpi_comm);
} else {
    MPI_Scatterv(mg_sdm[l].x, cnt_gatherv_max_ptr, disps_gatherv_max_ptr, *ddtype_gatherv_max_ptr, mg_sdm[l].x, 1, *ddtype_send_max_ptr, 0, comm_max_ptr->mpi_comm);
}


                // Loop for prolongation from aggregation level down to level 1
                for (int l = *lv_aggr_max_ptr-1; l >= *lv_aggr_med_ptr; l--) {
                    geometry_halocell_update_selectively(mg_sdm[l+1].x, &mg_sdm[l+1], mg_sdm[l+1].is_aggregated);

                    multigrid_prolongation_linear_on_nonuniform_grid(mg_sdm[l].r, mg_sdm[l+1].x,
                                                     &mg_sdm[l], &mg_sdm[l+1], l);
    

                    nx = mg_sdm[l].nx;
                    ny = mg_sdm[l].ny;
                    nz = mg_sdm[l].nz;
                    int size = (nx+2)*(ny+2)*(nz+2);
                    for (i = 0; i < size; i++) {
                        mg_sdm[l].x[i] += mg_sdm[l].r[i];
                    }

                    rbgs_iterator_poisson_matrix(mg_sdm[l].x, &mg_a_poisson[l], mg_sdm[l].b,
                                 &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                }
            }


l = *lv_aggr_med_ptr;
// #ifdef MPI_IN_PLACE
if (myrank == 0) {
    MPI_Scatterv(mg_sdm[l].x, cnt_gatherv_med_ptr, disps_gatherv_med_ptr, *ddtype_gatherv_med_ptr, MPI_IN_PLACE, 0, *ddtype_send_med_ptr, 0, comm_med_ptr->mpi_comm);
} else {
    MPI_Scatterv(mg_sdm[l].x, cnt_gatherv_med_ptr, disps_gatherv_med_ptr, *ddtype_gatherv_med_ptr, mg_sdm[l].x, 1, *ddtype_send_med_ptr, 0, comm_med_ptr->mpi_comm);
}


                // Loop for prolongation from aggregation level down to level 1
                for (int l = *lv_aggr_med_ptr-1; l >= 1; l--) {
                    geometry_halocell_update_selectively(mg_sdm[l+1].x, &mg_sdm[l+1], mg_sdm[l+1].is_aggregated);

                    multigrid_prolongation_linear_on_nonuniform_grid(mg_sdm[l].r, mg_sdm[l+1].x,
                                                     &mg_sdm[l], &mg_sdm[l+1], l);
    

                    nx = mg_sdm[l].nx;
                    ny = mg_sdm[l].ny;
                    nz = mg_sdm[l].nz;
                    int size = (nx+2)*(ny+2)*(nz+2);
                    for (i = 0; i < size; i++) {
                        mg_sdm[l].x[i] += mg_sdm[l].r[i];
                    }

                    rbgs_iterator_poisson_matrix(mg_sdm[l].x, &mg_a_poisson[l], mg_sdm[l].b,
                                 &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                }

                // Prolongation to actual solution
                geometry_halocell_update_selectively(mg_sdm[1].x, &mg_sdm[1], mg_sdm[1].is_aggregated);

                multigrid_prolongation_linear_on_nonuniform_grid(rsd, mg_sdm[1].x, sdm, &mg_sdm[1], 0);

                nx = sdm->nx;
                ny = sdm->ny;
                nz = sdm->nz;
                int size = (nx+2)*(ny+2)*(nz+2);
                for (i = 0; i < size; i++) {
                    sol[i] += rsd[i];
                }

                rbgs_iterator_poisson_matrix(sol, a_poisson, rhs, sdm, maxiteration, omega_sor, sdm->is_aggregated);

                multigrid_residual(rsd, a_poisson, sol, rhs, sdm, sdm->is_aggregated);

                vv_dot_3d_matrix(&rsd_val, rsd, rsd, nx, ny, nz, sdm->is_aggregated);

                if (myrank == 0) {
                    printf("[MG] cycle = %d, Error: %e %e %e\n", cyc, sqrt(rsd_val), sqrt(res0tol), sqrt(rsd_val/res0tol));
                }

                if (sqrt(rsd_val/res0tol) < tolerance) break;


            }

            if(myrank==0) printf("[MG] Total %d V-cycles end\n", cyc);

            free(rsd);
            rsd = NULL;


            break;

        case 3:
            // 3 aggregation level

            multigrid_residual(rsd, a_poisson, sol, rhs, sdm, sdm->is_aggregated);
            vv_dot_3d_matrix(&res0tol, rsd, rsd, nx, ny, nz, sdm->is_aggregated);

            for(cyc = 1; cyc <= n_vcycles; cyc++)
            {
                rbgs_iterator_poisson_matrix(sol, a_poisson, rhs, sdm, maxiteration, omega_sor, sdm->is_aggregated);
                multigrid_residual(rsd, a_poisson, sol, rhs, sdm, sdm->is_aggregated);
                multigrid_restriction(mg_sdm[1].b, rsd, &mg_sdm[1], sdm, 0);
                if(myrank==0) printf("[MG] Restriction from level 0 to level 1\n");

                for(l=1; l<=*lv_aggr_min_ptr-1; l++)
                {
                    size_t level_size = (mg_sdm[l].nx + 2) * (mg_sdm[l].ny + 2) * (mg_sdm[l].nz + 2);
                    memset(mg_sdm[l].x, 0, level_size * sizeof(double));
                    rbgs_iterator_poisson_matrix(mg_sdm[l].x, &mg_a_poisson[l], mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                    multigrid_residual(mg_sdm[l].r, &mg_a_poisson[l], mg_sdm[l].x, mg_sdm[l].b, &mg_sdm[l], mg_sdm[l].is_aggregated);
                    multigrid_restriction(mg_sdm[l+1].b, mg_sdm[l].r, &mg_sdm[l+1], &mg_sdm[l], l);
                    if(myrank==0) printf("[MG] Restriction from level %d to level %d\n", l, l+1);
                }

l = *lv_aggr_min_ptr;

if (myrank == 0) {
    MPI_Gatherv(MPI_IN_PLACE, 0, *ddtype_send_min_ptr, mg_sdm[l].b, cnt_gatherv_min_ptr, disps_gatherv_min_ptr, *ddtype_gatherv_min_ptr, 0, comm_min_ptr->mpi_comm);
} else {
    MPI_Gatherv(mg_sdm[l].b, 1, *ddtype_send_min_ptr, mg_sdm[l].b, cnt_gatherv_min_ptr, disps_gatherv_min_ptr, *ddtype_gatherv_min_ptr, 0, comm_min_ptr->mpi_comm);
}
                if ((*comm_min_ptr).myrank == 0) {

                    for (int l = *lv_aggr_min_ptr; l < *lv_aggr_med_ptr; l++) {
                        size_t level_size = (mg_sdm[l].nx + 2) * (mg_sdm[l].ny + 2) * (mg_sdm[l].nz + 2);
                        memset(mg_sdm[l].x, 0, level_size * sizeof(double));
                        rbgs_iterator_poisson_matrix(mg_sdm[l].x, &mg_a_poisson[l], mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                        multigrid_residual(mg_sdm[l].r, &mg_a_poisson[l], mg_sdm[l].x, mg_sdm[l].b, &mg_sdm[l], mg_sdm[l].is_aggregated);
                        multigrid_restriction(mg_sdm[l+1].b, mg_sdm[l].r, &mg_sdm[l+1], &mg_sdm[l], l);
                        printf("[MG] Restriction from level %d to level %d\n", l, l+1);
                    }
    
l = *lv_aggr_med_ptr;

if (myrank == 0) {
    MPI_Gatherv(MPI_IN_PLACE, 0, *ddtype_send_med_ptr, mg_sdm[l].b, cnt_gatherv_med_ptr, disps_gatherv_med_ptr, *ddtype_gatherv_med_ptr, 0, comm_med_ptr->mpi_comm);
} else {
    MPI_Gatherv(mg_sdm[l].b, 1, *ddtype_send_med_ptr, mg_sdm[l].b, cnt_gatherv_med_ptr, disps_gatherv_med_ptr, *ddtype_gatherv_med_ptr, 0, comm_med_ptr->mpi_comm);
}

                if ((*comm_med_ptr).myrank == 0) {
                    // Restriction phase
                    for (int l = *lv_aggr_med_ptr; l < *lv_aggr_max_ptr; l++) {
                        size_t level_size = (mg_sdm[l].nx + 2) * (mg_sdm[l].ny + 2) * (mg_sdm[l].nz + 2);
                        memset(mg_sdm[l].x, 0, level_size * sizeof(double));
                        rbgs_iterator_poisson_matrix(mg_sdm[l].x, &mg_a_poisson[l], mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                        multigrid_residual(mg_sdm[l].r, &mg_a_poisson[l], mg_sdm[l].x, mg_sdm[l].b, &mg_sdm[l], mg_sdm[l].is_aggregated);
                        multigrid_restriction(mg_sdm[l+1].b, mg_sdm[l].r, &mg_sdm[l+1], &mg_sdm[l], l);
                        printf("[MG] Restriction from level %d to level %d\n", l, l+1);
                    }

l = *lv_aggr_max_ptr;
  
if (myrank == 0) {
    MPI_Gatherv(MPI_IN_PLACE, 0, *ddtype_send_max_ptr, mg_sdm[l].b, cnt_gatherv_max_ptr, disps_gatherv_max_ptr, *ddtype_gatherv_max_ptr, 0, comm_max_ptr->mpi_comm);
} else {
    MPI_Gatherv(mg_sdm[l].b, 1, *ddtype_send_max_ptr, mg_sdm[l].b, cnt_gatherv_max_ptr, disps_gatherv_max_ptr, *ddtype_gatherv_max_ptr, 0, comm_max_ptr->mpi_comm);
}

                if ((*comm_max_ptr).myrank == 0) {
                    // Restriction phase
                    for (int l = *lv_aggr_max_ptr; l < n_levels; l++) {
                        size_t level_size = (mg_sdm[l].nx + 2) * (mg_sdm[l].ny + 2) * (mg_sdm[l].nz + 2);
                        memset(mg_sdm[l].x, 0, level_size * sizeof(double));
                        rbgs_iterator_poisson_matrix(mg_sdm[l].x, &mg_a_poisson[l], mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                        multigrid_residual(mg_sdm[l].r, &mg_a_poisson[l], mg_sdm[l].x, mg_sdm[l].b, &mg_sdm[l], mg_sdm[l].is_aggregated);
                        multigrid_restriction(mg_sdm[l+1].b, mg_sdm[l].r, &mg_sdm[l+1], &mg_sdm[l], l);
                        printf("[MG] Restriction from level %d to level %d\n", l, l+1);
                    }

                    // Solve at coarsest level
                    multigrid_solve_coarset_level(mg_sdm[n_levels].x, &mg_a_poisson[n_levels], mg_sdm[n_levels].b, &mg_sdm[n_levels], 1000, tolerance, omega_sor, mg_sdm[n_levels].is_aggregated);

                    multigrid_residual(mg_sdm[n_levels].r, &mg_a_poisson[n_levels], mg_sdm[n_levels].x, mg_sdm[n_levels].b, &mg_sdm[n_levels], mg_sdm[n_levels].is_aggregated);

                    printf("[MG] Solution in the coarsest level : x(1,1,1) = %.10e, residue = %.10e\n", mg_sdm[n_levels].x[IDX(1,1,1,mg_sdm[n_levels].nx,mg_sdm[n_levels].ny,mg_sdm[n_levels].nz)], mg_sdm[n_levels].r[IDX(1,1,1,mg_sdm[n_levels].nx,mg_sdm[n_levels].ny,mg_sdm[n_levels].nz)]);

#ifdef DEBUG_COARSEST
    multigrid_common_print_coarsest_level_solution(cyc, &mg_sdm[n_levels]);
#endif

                    // Prolongation phase
                    for (l = n_levels-1; l >= *lv_aggr_max_ptr; l--) {
                        multigrid_prolongation_linear_on_nonuniform_grid(mg_sdm[l].r, mg_sdm[l+1].x, &mg_sdm[l], &mg_sdm[l+1], l);

                        // x = x + r
                        nx = mg_sdm[l].nx;
                        ny = mg_sdm[l].ny;
                        nz = mg_sdm[l].nz;
                        for (i = 0; i < (nx+2)*(ny+2)*(nz+2); i++) {
                            mg_sdm[l].x[i] += mg_sdm[l].r[i];
                        }

                    rbgs_iterator_poisson_matrix(mg_sdm[l].x, &mg_a_poisson[l], mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                    }
                }

l = *lv_aggr_max_ptr;
// #ifdef MPI_IN_PLACE
if (myrank == 0) {
    MPI_Scatterv(mg_sdm[l].x, cnt_gatherv_max_ptr, disps_gatherv_max_ptr, *ddtype_gatherv_max_ptr, MPI_IN_PLACE, 0, *ddtype_send_max_ptr, 0, comm_max_ptr->mpi_comm);
} else {
    MPI_Scatterv(mg_sdm[l].x, cnt_gatherv_max_ptr, disps_gatherv_max_ptr, *ddtype_gatherv_max_ptr, mg_sdm[l].x, 1, *ddtype_send_max_ptr, 0, comm_max_ptr->mpi_comm);
}


                // Loop for prolongation from aggregation level down to level 1
                for (int l = *lv_aggr_max_ptr-1; l >= *lv_aggr_med_ptr; l--) {
                    geometry_halocell_update_selectively(mg_sdm[l+1].x, &mg_sdm[l+1], mg_sdm[l+1].is_aggregated);

                    multigrid_prolongation_linear_on_nonuniform_grid(mg_sdm[l].r, mg_sdm[l+1].x,
                                                     &mg_sdm[l], &mg_sdm[l+1], l);
    

                    nx = mg_sdm[l].nx;
                    ny = mg_sdm[l].ny;
                    nz = mg_sdm[l].nz;
                    int size = (nx+2)*(ny+2)*(nz+2);
                    for (i = 0; i < size; i++) {
                        mg_sdm[l].x[i] += mg_sdm[l].r[i];
                    }

                    rbgs_iterator_poisson_matrix(mg_sdm[l].x, &mg_a_poisson[l], mg_sdm[l].b,
                                 &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                }
            }


l = *lv_aggr_med_ptr;
// #ifdef MPI_IN_PLACE
if (myrank == 0) {
    MPI_Scatterv(mg_sdm[l].x, cnt_gatherv_med_ptr, disps_gatherv_med_ptr, *ddtype_gatherv_med_ptr, MPI_IN_PLACE, 0, *ddtype_send_med_ptr, 0, comm_med_ptr->mpi_comm);
} else {
    MPI_Scatterv(mg_sdm[l].x, cnt_gatherv_med_ptr, disps_gatherv_med_ptr, *ddtype_gatherv_med_ptr, mg_sdm[l].x, 1, *ddtype_send_med_ptr, 0, comm_med_ptr->mpi_comm);
}
                for (int l = *lv_aggr_med_ptr-1; l >= *lv_aggr_min_ptr; l--) {
                    geometry_halocell_update_selectively(mg_sdm[l+1].x, &mg_sdm[l+1], mg_sdm[l+1].is_aggregated);

                    multigrid_prolongation_linear_on_nonuniform_grid(mg_sdm[l].r, mg_sdm[l+1].x,
                                                     &mg_sdm[l], &mg_sdm[l+1], l);
    

                    nx = mg_sdm[l].nx;
                    ny = mg_sdm[l].ny;
                    nz = mg_sdm[l].nz;
                    int size = (nx+2)*(ny+2)*(nz+2);
                    for (i = 0; i < size; i++) {
                        mg_sdm[l].x[i] += mg_sdm[l].r[i];
                    }

                    rbgs_iterator_poisson_matrix(mg_sdm[l].x, &mg_a_poisson[l], mg_sdm[l].b,
                                 &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                }
            }

l = *lv_aggr_min_ptr;
// #ifdef MPI_IN_PLACE
if (myrank == 0) {
    MPI_Scatterv(mg_sdm[l].x, cnt_gatherv_min_ptr, disps_gatherv_min_ptr, *ddtype_gatherv_min_ptr, MPI_IN_PLACE, 0, *ddtype_send_min_ptr, 0, comm_min_ptr->mpi_comm);
} else {
    MPI_Scatterv(mg_sdm[l].x, cnt_gatherv_min_ptr, disps_gatherv_min_ptr, *ddtype_gatherv_min_ptr, mg_sdm[l].x, 1, *ddtype_send_min_ptr, 0, comm_min_ptr->mpi_comm);
}

                // Loop for prolongation from aggregation level down to level 1
                for (int l = *lv_aggr_min_ptr-1; l >= 1; l--) {
                    geometry_halocell_update_selectively(mg_sdm[l+1].x, &mg_sdm[l+1], mg_sdm[l+1].is_aggregated);

                    multigrid_prolongation_linear_on_nonuniform_grid(mg_sdm[l].r, mg_sdm[l+1].x,
                                                     &mg_sdm[l], &mg_sdm[l+1], l);
    

                    nx = mg_sdm[l].nx;
                    ny = mg_sdm[l].ny;
                    nz = mg_sdm[l].nz;
                    int size = (nx+2)*(ny+2)*(nz+2);
                    for (i = 0; i < size; i++) {
                        mg_sdm[l].x[i] += mg_sdm[l].r[i];
                    }

                    rbgs_iterator_poisson_matrix(mg_sdm[l].x, &mg_a_poisson[l], mg_sdm[l].b,
                                 &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                }

                // Prolongation to actual solution
                geometry_halocell_update_selectively(mg_sdm[1].x, &mg_sdm[1], mg_sdm[1].is_aggregated);

                multigrid_prolongation_linear_on_nonuniform_grid(rsd, mg_sdm[1].x, sdm, &mg_sdm[1], 0);

                nx = sdm->nx;
                ny = sdm->ny;
                nz = sdm->nz;
                int size = (nx+2)*(ny+2)*(nz+2);
                for (i = 0; i < size; i++) {
                    sol[i] += rsd[i];
                }

                rbgs_iterator_poisson_matrix(sol, a_poisson, rhs, sdm, maxiteration, omega_sor, sdm->is_aggregated);

                multigrid_residual(rsd, a_poisson, sol, rhs, sdm, sdm->is_aggregated);

                vv_dot_3d_matrix(&rsd_val, rsd, rsd, nx, ny, nz, sdm->is_aggregated);

                if (myrank == 0) {
                    printf("[MG] cycle = %d, Error: %e %e %e\n", cyc, sqrt(rsd_val), sqrt(res0tol), sqrt(rsd_val/res0tol));
                }

                if (sqrt(rsd_val/res0tol) < tolerance) break;


            }

            if(myrank==0) printf("[MG] Total %d V-cycles end\n", cyc);

            free(rsd);
            rsd = NULL;


            break;

        default:
            if (myrank == 0) printf("[Error] level_case should be 1, 2, or 3. Current: %d\n", level_case);
            MPI_Finalize();
            exit(EXIT_FAILURE);
    }


    

}

