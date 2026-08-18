#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>

#include "multigrid.h"
#include "poisson_matrix_operator.h" 
// #include "timer.h" 

#ifndef RESTRICTION_BLOCK_1D
#define RESTRICTION_BLOCK_1D 256
#endif

#ifndef RESTRICTION_BLOCK_X
#define RESTRICTION_BLOCK_X 8
#define RESTRICTION_BLOCK_Y 8
#define RESTRICTION_BLOCK_Z 4
#endif

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

#ifndef IDX
#define IDX(i,j,k,ni,nj) ((i)*(ni) + (j)*(nj) + (k))
#endif

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#define CUDA_CHECK(call) do {                                      \
    cudaError_t err__ = (call);                                    \
    if (err__ != cudaSuccess) {                                    \
        fprintf(stderr, "CUDA error %s:%d: %s\n",                  \
                __FILE__, __LINE__, cudaGetErrorString(err__));    \
        exit(EXIT_FAILURE);                                        \
    }                                                             \
} while (0)

/* Reusable work arrays avoid cudaMalloc/cudaFree in every V-cycle. */
static double *mg_restrict_dxf = NULL, *mg_restrict_dyf = NULL, *mg_restrict_dzf = NULL;
static size_t mg_restrict_dxf_capacity = 0, mg_restrict_dyf_capacity = 0, mg_restrict_dzf_capacity = 0;
static double *mg_prolong_dxp = NULL, *mg_prolong_dxn = NULL;
static double *mg_prolong_dyp = NULL, *mg_prolong_dyn = NULL;
static double *mg_prolong_dzp = NULL, *mg_prolong_dzn = NULL;
static size_t mg_prolong_x_capacity = 0, mg_prolong_y_capacity = 0, mg_prolong_z_capacity = 0;

static void multigrid_workspace_grow_gpu(double **buffer, size_t *capacity, size_t required)
{
    if (required <= *capacity) return;
    if (*buffer) CUDA_CHECK(cudaFree(*buffer));
    CUDA_CHECK(cudaMalloc((void **)buffer, required * sizeof(double)));
    *capacity = required;
}

static void multigrid_workspace_release_gpu(void)
{
    if (mg_restrict_dxf) CUDA_CHECK(cudaFree(mg_restrict_dxf));
    if (mg_restrict_dyf) CUDA_CHECK(cudaFree(mg_restrict_dyf));
    if (mg_restrict_dzf) CUDA_CHECK(cudaFree(mg_restrict_dzf));
    if (mg_prolong_dxp) CUDA_CHECK(cudaFree(mg_prolong_dxp));
    if (mg_prolong_dxn) CUDA_CHECK(cudaFree(mg_prolong_dxn));
    if (mg_prolong_dyp) CUDA_CHECK(cudaFree(mg_prolong_dyp));
    if (mg_prolong_dyn) CUDA_CHECK(cudaFree(mg_prolong_dyn));
    if (mg_prolong_dzp) CUDA_CHECK(cudaFree(mg_prolong_dzp));
    if (mg_prolong_dzn) CUDA_CHECK(cudaFree(mg_prolong_dzn));
    mg_restrict_dxf = mg_restrict_dyf = mg_restrict_dzf = NULL;
    mg_prolong_dxp = mg_prolong_dxn = NULL;
    mg_prolong_dyp = mg_prolong_dyn = NULL;
    mg_prolong_dzp = mg_prolong_dzn = NULL;
    mg_restrict_dxf_capacity = mg_restrict_dyf_capacity = mg_restrict_dzf_capacity = 0;
    mg_prolong_x_capacity = mg_prolong_y_capacity = mg_prolong_z_capacity = 0;
}


extern int lv_gdm_coarsest_x;
extern int lv_gdm_coarsest_y;
extern int lv_gdm_coarsest_z;

void multigrid_create_gpu(subdomain *sdm, int nlevel, int ncycle, int aggr_method, int aggr_level)
{
    int l;
    
    n_levels = nlevel;
    n_vcycles = ncycle;
    lv_aggregation = aggr_level;
    aggregation_type = aggr_method;

    if (n_levels <= 1) {
        // if (myrank == 0) {
        //     printf("[Error] The number of levels should be larger than 1. Current number: %d\n",
        //            n_levels);
        // }
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    mg_sdm = (subdomain *)malloc((n_levels + 1) * sizeof(subdomain));
    if (!mg_sdm) {
       if (myrank == 0) {
            printf("[Error] Failed to allocate mg_sdm\n");
        }
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    switch (aggregation_type)
    {
        case 0:
            multigrid_subdomain_create_no_aggregation_gpu(sdm);
            break;

        case 1:
            multigrid_subdomain_create_single_aggregation_gpu(sdm); 
            break;

        case 2:
            multigrid_subdomain_create_adaptive_aggregation_gpu(sdm);
            break;

        default:
            if (myrank == 0)
            {
                printf("[Error] Aggregation method should be 0, 1, or 2. Current: %d\n", aggregation_type);
            }
            MPI_Finalize();
            exit(EXIT_FAILURE);
    }

    multigrid_allocate_subdomain_variables_gpu();

    // if (myrank == 0) {
    //     printf("[MG-GPU] Aggregation info. of subdomain: %c %c %c\n",
    //        sdm->is_aggregated[0] ? 'T' : 'F',
    //        sdm->is_aggregated[1] ? 'T' : 'F',
    //        sdm->is_aggregated[2] ? 'T' : 'F');
    // }

    // for (int l = 1; l <= n_levels; l++) {
    //     if (myrank == 0) {
    //         printf("[MG] Aggregation info. of level %d: %c %c %c\n",
    //                l,
    //                mg_sdm[l].is_aggregated[0] ? 'T' : 'F',
    //                mg_sdm[l].is_aggregated[1] ? 'T' : 'F',
    //                mg_sdm[l].is_aggregated[2] ? 'T' : 'F');
    //     }
    // }

    // if (myrank == 0) {
    //     printf("[MG] Multigrid geometry constructed.\n");
    // }

    mg_a_poisson = (matrix_poisson *)malloc((n_levels + 1) * sizeof(matrix_poisson));
    if (!mg_a_poisson) {
        if (myrank == 0) {
            printf("[Error] Failed to allocate mg_a_poisson\n");
        }
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    for (int l = 1; l <= n_levels; l++) {
        matrix_poisson_create_gpu(&mg_a_poisson[l], &mg_sdm[l]);
    }

    // if (myrank == 0) {
    //     printf("[MG] Poisson matrix in multigrid constructed.\n");
    // }
}

void multigrid_subdomain_create_no_aggregation_gpu(subdomain *sdm)
{
    int l;
    int nx, ny, nz;
    double *dxm_f, *dxg_f, *xg_f;
    double *dym_f, *dyg_f, *yg_f;
    double *dzm_f, *dzg_f, *zg_f;

    // if (myrank == 0) {
    //     printf("[MG] Grid coarsening without aggregation.\n");
    // }

    if (lv_aggregation != 0) {
        if (myrank == 0) {
            printf("[Error] Aggregation level should be 0 for no aggregation method.\n");
        }
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    nx = sdm->nx;
    ny = sdm->ny;
    nz = sdm->nz;

    for (l = 1; l <= n_levels; l++) {
        mg_sdm[l].is_aggregated[0] = (comm_1d_x.nprocs == 1);
        mg_sdm[l].is_aggregated[1] = (comm_1d_y.nprocs == 1);
        mg_sdm[l].is_aggregated[2] = (comm_1d_z.nprocs == 1);
    }

    for (l = 1; l <= n_levels; l++) {
        // if (myrank == 0) printf("[MG] Grid coarsening at level %d\n", l);

        if (nx % 2 == 0) {
            int nx_old = nx;
            nx = nx / 2;
            lv_gdm_coarsest_x = l;
            // if (myrank == 0)
            //     printf("[MG] X-grid coarsening at level %d: %d grids reduced to %d grids.\n",
            //            l, nx_old, nx);
        } else {
            if (mg_sdm[l].is_aggregated[0]) {
                // if (myrank == 0) {
                //     printf("[MG] No more X-grid coarsening from level %d. Keeping X-grid number.\n", l);
                //     printf("[MG] The number of grid is %d and number processes is %d\n",
                //            nx, comm_1d_x.nprocs);
                // }
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

        if (ny % 2 == 0) {
            int ny_old = ny;
            ny = ny / 2;
            lv_gdm_coarsest_y = l;
            // if (myrank == 0)
            //     printf("[MG] Y-grid coarsening at level %d: %d grids reduced to %d grids.\n",
            //            l, ny_old, ny);
        } else {
            if (mg_sdm[l].is_aggregated[1]) {
                // if (myrank == 0) {
                //     printf("[MG] No more Y-grid coarsening from level %d. Keeping Y-grid number.\n", l);
                //     printf("[MG] The number of grid is %d and number processes is %d\n",
                //            ny, comm_1d_y.nprocs);
                // }
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

        if (nz % 2 == 0) {
            int nz_old = nz;
            nz = nz / 2;
            lv_gdm_coarsest_z = l;
            // if (myrank == 0)
            //     printf("[MG] Z-grid coarsening at level %d: %d grids reduced to %d grids.\n", l, nz_old, nz);
        } else {
            if (mg_sdm[l].is_aggregated[2]) {
                // if (myrank == 0) {
                //     printf("[MG] No more Z-grid coarsening from level %d. Keeping Z-grid number.\n", l);
                //     printf("[MG] The number of grid is %d and number processes is %d\n",
                //            nz, comm_1d_z.nprocs);
                // }
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

    // if (myrank == 0) {
    //     printf("[MG] Number of levels : %d\n", n_levels);
    //     printf("[MG] Final coarsest levels max, x, y, z directions : %d %d %d %d\n",
    //            lv_gdm_coarsest_max,
    //            lv_gdm_coarsest_x,
    //            lv_gdm_coarsest_y,
    //            lv_gdm_coarsest_z);
    // }

    if (n_levels > lv_gdm_coarsest_max) {
        // if (myrank == 0) {
        //     printf("[MG] Number of levels is greater than the coarsest level in global domain.\n");
        //     printf("[MG] It is not allowed and reduce the number of levels.\n");
        // }
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    // if (myrank == 0) printf("[MG] Generating grid dimension\n");

    dxm_f = sdm->dxm;
    dxg_f = sdm->dxg;
    xg_f  = sdm->xg;

    for (l = 1; l <= n_levels; l++) {
        nx = mg_sdm[l].nx;

        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].dxm, (nx + 2) * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].dxg, (nx + 2) * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].xg,  (nx + 2) * sizeof(double)));

        multigrid_subdomain_no_aggregation_make_grid_gpu(
            mg_sdm[l].dxm, mg_sdm[l].dxg, mg_sdm[l].xg,
            nx,
            dxm_f, dxg_f, xg_f,
            sdm->ox,
            l,
            lv_gdm_coarsest_x,
            comm_1d_x,
            'x');

        dxm_f = mg_sdm[l].dxm;
        dxg_f = mg_sdm[l].dxg;
        xg_f  = mg_sdm[l].xg;
    }

    dym_f = sdm->dym;
    dyg_f = sdm->dyg;
    yg_f  = sdm->yg;

    for (l = 1; l <= n_levels; l++) {
        ny = mg_sdm[l].ny;

        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].dym, (ny + 2) * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].dyg, (ny + 2) * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].yg,  (ny + 2) * sizeof(double)));

        multigrid_subdomain_no_aggregation_make_grid_gpu(
            mg_sdm[l].dym, mg_sdm[l].dyg, mg_sdm[l].yg,
            ny,
            dym_f, dyg_f, yg_f,
            sdm->oy,
            l,
            lv_gdm_coarsest_y,
            comm_1d_y,
            'y');

        dym_f = mg_sdm[l].dym;
        dyg_f = mg_sdm[l].dyg;
        yg_f  = mg_sdm[l].yg;
    }

    dzm_f = sdm->dzm;
    dzg_f = sdm->dzg;
    zg_f  = sdm->zg;

    for (l = 1; l <= n_levels; l++) {
        nz = mg_sdm[l].nz;

        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].dzm, (nz + 2) * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].dzg, (nz + 2) * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].zg,  (nz + 2) * sizeof(double)));

        multigrid_subdomain_no_aggregation_make_grid_gpu(
            mg_sdm[l].dzm, mg_sdm[l].dzg, mg_sdm[l].zg,
            nz, dzm_f, dzg_f, zg_f,
            sdm->oz,
            l, lv_gdm_coarsest_z, comm_1d_z, 'z');
            
        dzm_f = mg_sdm[l].dzm;
        dzg_f = mg_sdm[l].dzg;
        zg_f  = mg_sdm[l].zg;
    }
}


void multigrid_subdomain_create_single_aggregation_gpu(subdomain *sdm) 
{
    int l;
    int nx, ny, nz;
    int nx_lv_aggregation = 0, ny_lv_aggregation = 0, nz_lv_aggregation = 0;

    double *dxm_f, *dxg_f, *xg_f;
    double *dym_f, *dyg_f, *yg_f;
    double *dzm_f, *dzg_f, *zg_f;

    // if (myrank == 0) {printf("[MG] Grid coarsening with aggretation level = %d\n", lv_aggregation);}

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
        // if (myrank == 0) printf("[MG] Grid coarsening at level %d\n", l);

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
                // if (myrank == 0)
                // {
                //     printf("[MG] X-grid coarsening at level %d: %d grids reduced to %d grids.\n", l, nx, nx / 2);
                // }
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
                // if (myrank == 0)
                // {
                //     printf("[MG] Y-grid coarsening at level %d: %d grids reduced to %d grids.\n", l, ny, ny / 2);
                // }
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
                // if (myrank == 0)
                // {
                //     printf("[MG] Z-grid coarsening at level %d: %d grids reduced to %d grids.\n", l, nz, nz / 2);
                // }
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
                // if (myrank == 0)
                // {
                //     printf("[MG] X-grid coarsening at level %d: %d grids reduced to %d grids.\n", l, nx, nx / 2);
                //     printf("[MG] X-grid aggregation at level %d: %d grids aggregated to %d grids.\n", l, nx / 2, mg_sdm[l].nx);
                // }
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
                // if (myrank == 0) {
                //     printf("[MG] Y-grid coarsening at level %d: %d grids reduced to %d grids.\n", l, ny, ny / 2);
                //     printf("[MG] Y-grid aggregation at level %d: %d grids aggregated to %d grids.\n", l, ny / 2, mg_sdm[l].ny);
                // }
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
                // if (myrank == 0) {
                //     printf("[MG] Z-grid coarsening at level %d: %d grids reduced to %d grids.\n", l, nz, nz / 2);
                //     printf("[MG] Z-grid aggregation at level %d: %d grids aggregated to %d grids.\n", l, nz / 2, mg_sdm[l].nz);
                // }
            }
        }
        else if (l > lv_aggregation)
        {
            // X
            if (nx % 2 == 1) 
            {
                // if (myrank == 0) 
                // {
                //     printf("[MG] No more X-grid coarsening from level %d. Keeping X-grid number.\n", l);
                //     printf("[MG] The number of grid is %d and number processes is %d\n", nx, comm_1d_x.nprocs);
                // }
                mg_sdm[l].nx = nx;
            } 
            else {
                // if (myrank == 0){
                //     printf("[MG] X-grid coarsening at level %d: %d grids reduced to %d grids.\n", l, nx, nx / 2);
                // }
                mg_sdm[l].nx = nx / 2;
                lv_gdm_coarsest_x = l;
            }
            // Y
            if (ny % 2 == 1) 
            {
                // if (myrank == 0) 
                // {
                //     printf("[MG] No more Y-grid coarsening from level %d. Keeping Y-grid number.\n", l);
                //     printf("[MG] The number of grid is %d and number processes is %d\n", ny, comm_1d_y.nprocs);
                // }
                mg_sdm[l].ny = ny;
            } 
            else {
                // if (myrank == 0){
                //     printf("[MG] Y-grid coarsening at level %d: %d grids reduced to %d grids.\n", l, ny, ny / 2);
                // }
                mg_sdm[l].ny = ny / 2;
                lv_gdm_coarsest_y = l;
            }
            // Z
            if (nz % 2 == 1) {
                // if (myrank == 0) 
                // {
                //     printf("[MG] No more Z-grid coarsening from level %d. Keeping Z-grid number.\n", l);
                //     printf("[MG] The number of grid is %d and number processes is %d\n", nz, comm_1d_z.nprocs);
                // }
                 mg_sdm[l].nz = nz;
            } 
            else {
                // if (myrank == 0){
                //     printf("[MG] Z-grid coarsening at level %d: %d grids reduced to %d grids.\n", l, nz, nz / 2);
                // }
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
    // if (myrank == 0) {
    //     printf("[MG] Number of levels : %d\n", n_levels);
    //     printf("[MG] Final coarsest levels max, x, y, z directions : %d %d %d %d\n",
    //            lv_gdm_coarsest_max, lv_gdm_coarsest_x, lv_gdm_coarsest_y, lv_gdm_coarsest_z);
    // }

    // Check the number of levels and the max coarsest level
    if (n_levels > lv_gdm_coarsest_max) {
        if (myrank == 0) {
            printf("[MG] Number of levels is greater than the coarsest level in global domain.\n");
            printf("[MG] It is not allowed and reduce the number of levels.\n");
        }
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    // if (myrank == 0){
    //     printf("[MG] Final gather level = %d\n",lv_aggregation);
    // }
    if (lv_aggregation > n_levels)
    {
        if (myrank == 0){
            printf("[MG] Gather level is greater than the coarsest level and no gather level exists\n");
        }
    }
    // if (myrank == 0) printf("[MG] Generating grid dimension\n");

    // In x-direction
    dxm_f = sdm->dxm;
    dxg_f = sdm->dxg;
    xg_f  = sdm->xg;

    for (l = 1; l <= n_levels; l++) {
        nx = mg_sdm[l].nx;
        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].dxm, (nx + 2) * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].dxg, (nx + 2) * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].xg, (nx + 2) * sizeof(double)));

        multigrid_subdomain_aggregation_make_grid_gpu(mg_sdm[l].dxm, mg_sdm[l].dxg, mg_sdm[l].xg, nx, dxm_f, dxg_f, xg_f, sdm->ox, l, lv_gdm_coarsest_x, lv_aggregation_x, nx_lv_aggregation, comm_1d_x, 'x');

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
        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].dym, (ny + 2) * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].dyg, (ny + 2) * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].yg, (ny + 2) * sizeof(double)));

        multigrid_subdomain_aggregation_make_grid_gpu(mg_sdm[l].dym, mg_sdm[l].dyg, mg_sdm[l].yg, ny, dym_f, dyg_f, yg_f, sdm->oy, l, lv_gdm_coarsest_y, lv_aggregation_y, ny_lv_aggregation, comm_1d_y, 'y');

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
        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].dzm, (nz + 2) * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].dzg, (nz + 2) * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].zg, (nz + 2) * sizeof(double)));

        multigrid_subdomain_aggregation_make_grid_gpu(mg_sdm[l].dzm, mg_sdm[l].dzg, mg_sdm[l].zg, nz, dzm_f, dzg_f, zg_f, sdm->oz, l, lv_gdm_coarsest_z, lv_aggregation_z, nz_lv_aggregation, comm_1d_z, 'z');

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

void multigrid_subdomain_create_adaptive_aggregation_gpu(subdomain *sdm) 
{
    int l;
    int nx, ny, nz;
    int nx_lv_aggregation = 0, ny_lv_aggregation = 0, nz_lv_aggregation = 0;

    double *dxm_f, *dxg_f, *xg_f;
    double *dym_f, *dyg_f, *yg_f;
    double *dzm_f, *dzg_f, *zg_f;

    // if (myrank == 0) {printf("[MG] Grid coarsening with adaptive aggretation. lv_aggregation is neglected\n");}

    nx = sdm->nx;
    ny = sdm->ny;
    nz = sdm->nz;

    lv_aggregation_x = n_levels+1;
    lv_aggregation_y = n_levels+1;
    lv_aggregation_z = n_levels+1;


    for (l = 1; l <= n_levels; l++) {
        // if (myrank == 0) printf("[MG] Grid coarsening at level %d\n", l);
    
        // X
        if (nx % 2 == 1) 
        {
            mg_sdm[l].nx = nx;
            // if (myrank == 0) 
            // {
            //     printf("[MG] No more x-grid coarsening at level %d in serial. Keeping x-grid number\n", l);
            //     printf("[MG] Aggregation level in x-direction is %d\n", lv_aggregation_x);
            //     printf("[MG] The number of grid is %d and number processes is %d\n", mg_sdm[l].nx, comm_1d_x.nprocs);
            // }
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
                // if (myrank == 0)
                // {
                //     printf("[MG] x-grid is odd number at level %d in parallel.\n", l);
                //     printf("[MG] Aggregation level in x-direction prescribed as %d\n", lv_aggregation_x);
                //     printf("[MG] The number of grid is %d and number processes is %d\n", mg_sdm[l].nx, comm_1d_x.nprocs);
                // }
            }
            else{
                mg_sdm[l].nx = nx/2;
                lv_gdm_coarsest_x = l;

                // if (myrank == 0) 
                // {
                //     printf("[MG] x-grid coarsening at level %d: %d grids reduced to %d grids.\n", l, nx, mg_sdm[l].nx);
                // }
            }  
        }

        // Y
        if (ny % 2 == 1) 
        {
            mg_sdm[l].ny = ny;
            // if (myrank == 0) 
            // {
            //     printf("[MG] No more y-grid coarsening at level %d in serial. Keeping y-grid number\n", l);
            //     printf("[MG] Aggregation level in y-direction is %d\n", lv_aggregation_y);
            //     printf("[MG] The number of grid is %d and number processes is %d\n", mg_sdm[l].ny, comm_1d_y.nprocs);
            // }
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
                // if (myrank == 0)
                // {
                //     printf("[MG] y-grid is odd number at level %d in parallel.\n", l);
                //     printf("[MG] Aggregation level in y-direction prescribed as %d\n", lv_aggregation_y);
                //     printf("[MG] The number of grid is %d and number processes is %d\n", mg_sdm[l].ny, comm_1d_y.nprocs);
                // }
            }
            else{
                mg_sdm[l].ny = ny/2;
                lv_gdm_coarsest_y = l;
                // if (myrank == 0) 
                // {
                //     printf("[MG] y-grid coarsening at level %d: %d grids reduced to %d grids.\n", l, ny, mg_sdm[l].ny);
                // }
            }  
        }

        // Z
        if (nz % 2 == 1) 
        {
            mg_sdm[l].nz = nz;
            // if (myrank == 0) 
            // {
            //     printf("[MG] No more z-grid coarsening at level %d in serial. Keeping z-grid number\n", l);
            //     printf("[MG] Aggregation level in z-direction is %d\n", lv_aggregation_z);
            //     printf("[MG] The number of grid is %d and number processes is %d\n", mg_sdm[l].nz, comm_1d_z.nprocs);
            // }
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
                // if (myrank == 0)
                // {
                //     printf("[MG] z-grid is odd number at level %d in parallel.\n", l);
                //     printf("[MG] Aggregation level in z-direction prescribed as %d\n", lv_aggregation_z);
                //     printf("[MG] The number of grid is %d and number processes is %d\n", mg_sdm[l].nz, comm_1d_z.nprocs);
                // }
            }
            else{
                mg_sdm[l].nz = nz/2;
                lv_gdm_coarsest_z = l;
                // if (myrank == 0) 
                // {
                //     printf("[MG] z-grid coarsening at level %d: %d grids reduced to %d grids.\n", l, nz, mg_sdm[l].nz);
                // }
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
    // if (myrank == 0) {
    //     printf("[MG] Number of levels : %d\n", n_levels);
    //     printf("[MG] Final aggregation levels max, x, y, z directions : %d %d %d %d\n",
    //            lv_aggregation_max, lv_aggregation_x, lv_aggregation_y, lv_aggregation_z);
    //     printf("[MG] Final coarsest levels max, x, y, z directions : %d %d %d %d\n",
    //            lv_gdm_coarsest_max, lv_gdm_coarsest_x, lv_gdm_coarsest_y, lv_gdm_coarsest_z);
    // }

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
        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].dxm, (nx + 2) * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].dxg, (nx + 2) * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].xg, (nx + 2) * sizeof(double)));

        multigrid_subdomain_aggregation_make_grid_gpu(mg_sdm[l].dxm, mg_sdm[l].dxg, mg_sdm[l].xg, nx, dxm_f, dxg_f, xg_f, sdm->ox, l, lv_gdm_coarsest_x, lv_aggregation_x, nx_lv_aggregation, comm_1d_x, 'x');

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
        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].dym, (ny + 2) * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].dyg, (ny + 2) * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].yg, (ny + 2) * sizeof(double)));

        multigrid_subdomain_aggregation_make_grid_gpu(mg_sdm[l].dym, mg_sdm[l].dyg, mg_sdm[l].yg, ny, dym_f, dyg_f, yg_f, sdm->oy, l, lv_gdm_coarsest_y, lv_aggregation_y, ny_lv_aggregation, comm_1d_y, 'y');

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
        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].dzm, (nz + 2) * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].dzg, (nz + 2) * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].zg, (nz + 2) * sizeof(double)));

        multigrid_subdomain_aggregation_make_grid_gpu(mg_sdm[l].dzm, mg_sdm[l].dzg, mg_sdm[l].zg, nz, dzm_f, dzg_f, zg_f, sdm->oz, l, lv_gdm_coarsest_z, lv_aggregation_z, nz_lv_aggregation, comm_1d_z, 'z');

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





__global__
static void make_grid_zero_kernel(double *dxm, double *dxg, double *xg, int nx)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i <= nx + 1) {
        dxm[i] = 0.0;
        dxg[i] = 0.0;
        xg[i]  = 0.0;
    }
}

__global__
static void make_grid_dxm_coarsen_kernel(double *dxm, const double *dxm_f, int nx)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;

    if (i <= nx) {
        dxm[i] = dxm_f[2*i - 1] + dxm_f[2*i];

    }

    if (i == 1) {
        dxm[0] = dxm_f[0];
        dxm[nx + 1] = dxm_f[2*nx + 1];
    }
}

__global__
static void make_grid_copy_kernel(double *dxm, double *dxg, double *xg, const double *dxm_f, const double *dxg_f, const double *xg_f, int nx)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i <= nx + 1) {
        dxm[i] = dxm_f[i];
        dxg[i] = dxg_f[i];
        xg[i]  = xg_f[i];
    }
}

__global__
static void make_grid_dxg_xg_kernel(double *dxm, double *dxg, double *xg, int nx, double ox)
{
    if (blockIdx.x == 0 && threadIdx.x == 0) {
        xg[0] = ox - 0.5 * dxm[0];

        for (int i = 1; i <= nx; i++) {
            dxg[i] = 0.5 * (dxm[i] + dxm[i-1]);
            xg[i]  = xg[i-1] + dxg[i];
        }

        dxg[nx + 1] = 0.5 * (dxm[nx + 1] + dxm[nx]);
        xg[nx + 1]  = xg[nx] + dxg[nx + 1];
    }
}


__global__
static void copy_one_kernel(double *dst, int dst_idx, const double *src, int src_idx)
{
    dst[dst_idx] = src[src_idx];
}

void multigrid_subdomain_no_aggregation_make_grid_gpu( double *dxm, double *dxg, double *xg, int nx, const double *dxm_f, const double *dxg_f, const double *xg_f, double ox, int lv_cur, int lv_coarsest, cart_comm_1d comm_1d, char dir)
{
    MPI_Request request1[2], request2[2];

    // if (myrank == 0) {
    //     printf("[MG] level : %d, dir : %c\n", lv_cur, dir);
    // }

    int block1d = 256;
    int grid1d = (nx + 2 + block1d - 1) / block1d;

    make_grid_zero_kernel<<<grid1d, block1d>>>(dxm, dxg, xg, nx);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    if (lv_cur <= lv_coarsest) {
        make_grid_dxm_coarsen_kernel<<<grid1d, block1d>>>(dxm, dxm_f, nx);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());

        if (comm_1d.west_rank != MPI_PROC_NULL) {
            MPI_Isend(dxm + 1, 1, MPI_DOUBLE,
                      comm_1d.west_rank, 1,
                      comm_1d.mpi_comm, &request1[0]);

            MPI_Irecv(dxm + 0, 1, MPI_DOUBLE,
                      comm_1d.west_rank, 2,
                      comm_1d.mpi_comm, &request1[1]);
        } else {
            copy_one_kernel<<<1, 1>>>(dxm, 0, dxm_f, 0);
            CUDA_CHECK(cudaGetLastError());
        }

        if (comm_1d.east_rank != MPI_PROC_NULL) {
            MPI_Isend(dxm + nx, 1, MPI_DOUBLE,
                      comm_1d.east_rank, 2,
                      comm_1d.mpi_comm, &request2[0]);

            MPI_Irecv(dxm + nx + 1, 1, MPI_DOUBLE,
                      comm_1d.east_rank, 1,
                      comm_1d.mpi_comm, &request2[1]);
        } else {
            copy_one_kernel<<<1, 1>>>(dxm, nx + 1, dxm_f, 2 * nx + 1);
            CUDA_CHECK(cudaGetLastError());
        }

        if (comm_1d.west_rank != MPI_PROC_NULL) {
            MPI_Waitall(2, request1, MPI_STATUSES_IGNORE);
        }

        if (comm_1d.east_rank != MPI_PROC_NULL) {
            MPI_Waitall(2, request2, MPI_STATUSES_IGNORE);
        }

        CUDA_CHECK(cudaDeviceSynchronize());

        make_grid_dxg_xg_kernel<<<1, 1>>>(dxm, dxg, xg, nx, ox);
        CUDA_CHECK(cudaGetLastError());
    } else {
        make_grid_copy_kernel<<<grid1d, block1d>>>(
            dxm, dxg, xg,
            dxm_f, dxg_f, xg_f,
            nx
        );
        CUDA_CHECK(cudaGetLastError());
    }

    CUDA_CHECK(cudaDeviceSynchronize());
    
  
}

__global__
static void make_grid_set_x0_from_ox_kernel(double *xg, const double *dxm, double ox)
{
    if (blockIdx.x == 0 && threadIdx.x == 0)
        xg[0] = ox - 0.5 * dxm[0];
}

__global__
static void make_grid_set_x0_after_aggregation_kernel(double *xg, const double *dxm, const double *dxm_f, const double *xg_f)
{
    if (blockIdx.x == 0 && threadIdx.x == 0)
        xg[0] = (xg_f[0] + 0.5 * dxm_f[0]) - 0.5 * dxm[0];
}

/* xg is a recurrence, so this retains the CPU order and round-off behavior. */
__global__
static void make_grid_finish_coordinates_kernel(double *dxm, double *dxg,
                                                 double *xg, int nx)
{
    if (blockIdx.x == 0 && threadIdx.x == 0) {
        for (int i = 1; i <= nx; ++i) {
            dxg[i] = 0.5 * (dxm[i] + dxm[i - 1]);
            xg[i] = xg[i - 1] + dxg[i];
        }
        dxg[nx + 1] = 0.5 * (dxm[nx + 1] + dxm[nx]);
        xg[nx + 1] = xg[nx] + dxg[nx + 1];
    }
}

void multigrid_subdomain_aggregation_make_grid_gpu(double *dxm, double *dxg, double *xg, int nx, const double *dxm_f, const double *dxg_f, const double *xg_f, double ox, int lv_cur, int lv_coarsest, int lv_aggregation, int nx_lv_aggregation, cart_comm_1d comm_1d, char dir)
{
    MPI_Request request1[2], request2[2];
    const int block1d = 256;
    const int grid1d = (nx + 2 + block1d - 1) / block1d;

    // if (myrank == 0) {
    //     if (lv_cur == lv_aggregation)
    //         printf("[MG-GPU] level(aggregation) : %d, dir : %c\n",
    //                lv_cur, dir);
    //     else
    //         printf("[MG-GPU] level : %d, dir : %c\n", lv_cur, dir);
    // }

    make_grid_zero_kernel<<<grid1d, block1d>>>(dxm, dxg, xg, nx);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    if (lv_cur > lv_coarsest) {
        make_grid_copy_kernel<<<grid1d, block1d>>>(
            dxm, dxg, xg, dxm_f, dxg_f, xg_f, nx);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
        return;
    }

    if (lv_cur < lv_aggregation) {
        make_grid_dxm_coarsen_kernel<<<grid1d, block1d>>>(dxm, dxm_f, nx);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());

        if (comm_1d.west_rank != MPI_PROC_NULL) {
            MPI_Isend(dxm + 1, 1, MPI_DOUBLE, comm_1d.west_rank, 1,
                      comm_1d.mpi_comm, &request1[0]);
            MPI_Irecv(dxm, 1, MPI_DOUBLE, comm_1d.west_rank, 2,
                      comm_1d.mpi_comm, &request1[1]);
        } else {
            copy_one_kernel<<<1, 1>>>(dxm, 0, dxm_f, 0);
            CUDA_CHECK(cudaGetLastError());
        }

        if (comm_1d.east_rank != MPI_PROC_NULL) {
            MPI_Isend(dxm + nx, 1, MPI_DOUBLE, comm_1d.east_rank, 2,
                      comm_1d.mpi_comm, &request2[0]);
            MPI_Irecv(dxm + nx + 1, 1, MPI_DOUBLE, comm_1d.east_rank, 1,
                      comm_1d.mpi_comm, &request2[1]);
        } else {
            copy_one_kernel<<<1, 1>>>(dxm, nx + 1, dxm_f, 2 * nx + 1);
            CUDA_CHECK(cudaGetLastError());
        }

        if (comm_1d.west_rank != MPI_PROC_NULL)
            MPI_Waitall(2, request1, MPI_STATUSES_IGNORE);
        if (comm_1d.east_rank != MPI_PROC_NULL)
            MPI_Waitall(2, request2, MPI_STATUSES_IGNORE);

        make_grid_set_x0_from_ox_kernel<<<1, 1>>>(xg, dxm, ox);
        make_grid_finish_coordinates_kernel<<<1, 1>>>(dxm, dxg, xg, nx);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
        return;
    }

    if (lv_cur == lv_aggregation) {
        /* Local coarse result has nx_lv_aggregation values plus two ghosts. */
        double *dxm_local = NULL;
        const int local_grid =
            (nx_lv_aggregation + 2 + block1d - 1) / block1d;
        CUDA_CHECK(cudaMalloc((void **)&dxm_local,
                              (size_t)(nx_lv_aggregation + 2) *
                                  sizeof(double)));
        CUDA_CHECK(cudaMemset(dxm_local, 0,
                              (size_t)(nx_lv_aggregation + 2) *
                                  sizeof(double)));

        make_grid_dxm_coarsen_kernel<<<local_grid, block1d>>>(
            dxm_local, dxm_f, nx_lv_aggregation);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());


        MPI_Allgather(dxm_local + 1, nx_lv_aggregation, MPI_DOUBLE,
                      dxm + 1, nx_lv_aggregation, MPI_DOUBLE,
                      comm_1d.mpi_comm);

        /* Each rank stages its local boundary; Bcast selects the true owner. */
        copy_one_kernel<<<1, 1>>>(dxm, 0, dxm_f, 0);
        copy_one_kernel<<<1, 1>>>(dxm, nx + 1, dxm_f,
                                  2 * nx_lv_aggregation + 1);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());

        MPI_Bcast(dxm, 1, MPI_DOUBLE, 0, comm_1d.mpi_comm);
        MPI_Bcast(dxm + nx + 1, 1, MPI_DOUBLE, comm_1d.nprocs - 1,
                  comm_1d.mpi_comm);

        make_grid_set_x0_from_ox_kernel<<<1, 1>>>(xg, dxm, ox);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
        MPI_Bcast(xg, 1, MPI_DOUBLE, 0, comm_1d.mpi_comm);

        make_grid_finish_coordinates_kernel<<<1, 1>>>(dxm, dxg, xg, nx);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
        CUDA_CHECK(cudaFree(dxm_local));
        return;
    }

    /* lv_cur > lv_aggregation: the grid is already globally aggregated. */
    make_grid_dxm_coarsen_kernel<<<grid1d, block1d>>>(dxm, dxm_f, nx);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    make_grid_set_x0_after_aggregation_kernel<<<1, 1>>>(
        xg, dxm, dxm_f, xg_f);
    make_grid_finish_coordinates_kernel<<<1, 1>>>(dxm, dxg, xg, nx);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
}





void multigrid_allocate_subdomain_variables_gpu()
{
    for (int l = 1; l <= n_levels; l++) {
        int nx = mg_sdm[l].nx;
        int ny = mg_sdm[l].ny;
        int nz = mg_sdm[l].nz;

        size_t total_size = (size_t)(nx + 2) * (size_t)(ny + 2) * (size_t)(nz + 2);

        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].b, total_size * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].r, total_size * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void**)&mg_sdm[l].x, total_size * sizeof(double)));

        CUDA_CHECK(cudaMemset(mg_sdm[l].b, 0, total_size * sizeof(double)));
        CUDA_CHECK(cudaMemset(mg_sdm[l].r, 0, total_size * sizeof(double)));
        CUDA_CHECK(cudaMemset(mg_sdm[l].x,  0, total_size * sizeof(double)));
    }

   
    for (int l = 1; l <= n_levels; l++) {
        geometry_subdomain_ddt_create(&mg_sdm[l]);
    }
}






__global__
static void build_dxf_kernel(double *dxf, const double *dxm_c, const double *dxm_f, int nx_f, int level, int lv_gdm_coarsest_x)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (level >= lv_gdm_coarsest_x) {
        if (i <= nx_f + 1) {
            dxf[i] = dxm_f[i];
        }
    } else {
        if (i >= 1 && i <= nx_f) {
            int i_c = (i - 1) / 2 + 1;
            dxf[i] = dxm_c[i_c] - dxm_f[i];
        }
    }
}

__global__
static void build_dyf_kernel(double *dyf, const double *dym_c, const double *dym_f, int ny_f, int level, int lv_gdm_coarsest_y)
{
    int j = blockIdx.x * blockDim.x + threadIdx.x;

    if (level >= lv_gdm_coarsest_y) {
        if (j <= ny_f + 1) {
            dyf[j] = dym_f[j];
        }
    } else {
        if (j >= 1 && j <= ny_f) {
            int j_c = (j - 1) / 2 + 1;
            dyf[j] = dym_c[j_c] - dym_f[j];
        }
    }
}

__global__
static void build_dzf_kernel(double *dzf, const double *dzm_c, const double *dzm_f, int nz_f, int level, int lv_gdm_coarsest_z)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;

    if (level >= lv_gdm_coarsest_z) {
        if (k <= nz_f + 1) {
            dzf[k] = dzm_f[k];
        }
    } else {
        if (k >= 1 && k <= nz_f) {
            int k_c = (k - 1) / 2 + 1;
            dzf[k] = dzm_c[k_c] - dzm_f[k];
        }
    }
}


__global__
static void restriction_kernel(double *val_c, const double *val_f, const double *dxf, const double *dyf, const double *dzf,
                               int nx_c, int ny_c, int nz_c,
                               int nx_f, int ny_f, int nz_f,
                               int dm_c_ny, int dm_c_nz,
                               int i_gl_c, int j_gl_c, int k_gl_c,
                               int level, int lv_gdm_coarsest_x, int lv_gdm_coarsest_y, int lv_gdm_coarsest_z)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

    if (i > nx_c || j > ny_c || k > nz_c) return;

    int i_stride_f = (level >= lv_gdm_coarsest_x) ? 1 : 2;
    int j_stride_f = (level >= lv_gdm_coarsest_y) ? 1 : 2;
    int k_stride_f = (level >= lv_gdm_coarsest_z) ? 1 : 2;

    int i_offset_f = (level >= lv_gdm_coarsest_x) ? 0 : 1;
    int j_offset_f = (level >= lv_gdm_coarsest_y) ? 0 : 1;
    int k_offset_f = (level >= lv_gdm_coarsest_z) ? 0 : 1;


    int i_c = i + i_gl_c;
    int j_c = j + j_gl_c;
    int k_c = k + k_gl_c;

    int ip_f = i * i_stride_f;
    int jp_f = j * j_stride_f;
    int kp_f = k * k_stride_f;

    int iz_f = ip_f - i_offset_f;
    int jz_f = jp_f - j_offset_f;
    int kz_f = kp_f - k_offset_f;

    int ni_c = (dm_c_ny + 2) * (dm_c_nz + 2);
    int nj_c = (dm_c_nz + 2);

    int ni_f = (ny_f + 2) * (nz_f + 2);
    int nj_f = (nz_f + 2);

    double v000 = dxf[iz_f] * dyf[jz_f] * dzf[kz_f];
    double v100 = dxf[ip_f] * dyf[jz_f] * dzf[kz_f];
    double v010 = dxf[iz_f] * dyf[jp_f] * dzf[kz_f];
    double v110 = dxf[ip_f] * dyf[jp_f] * dzf[kz_f];
    double v001 = dxf[iz_f] * dyf[jz_f] * dzf[kp_f];
    double v101 = dxf[ip_f] * dyf[jz_f] * dzf[kp_f];
    double v011 = dxf[iz_f] * dyf[jp_f] * dzf[kp_f];
    double v111 = dxf[ip_f] * dyf[jp_f] * dzf[kp_f];

    double vol_c = v000 + v100 + v010 + v110 + v001 + v101 + v011 + v111;

    val_c[IDX(i_c,j_c,k_c,ni_c,nj_c)] =
        ( v000 * val_f[IDX(ip_f,jp_f,kp_f,ni_f,nj_f)]
        + v100 * val_f[IDX(iz_f,jp_f,kp_f,ni_f,nj_f)]
        + v010 * val_f[IDX(ip_f,jz_f,kp_f,ni_f,nj_f)]
        + v110 * val_f[IDX(iz_f,jz_f,kp_f,ni_f,nj_f)]
        + v001 * val_f[IDX(ip_f,jp_f,kz_f,ni_f,nj_f)]
        + v101 * val_f[IDX(iz_f,jp_f,kz_f,ni_f,nj_f)]
        + v011 * val_f[IDX(ip_f,jz_f,kz_f,ni_f,nj_f)]
        + v111 * val_f[IDX(iz_f,jz_f,kz_f,ni_f,nj_f)] ) / vol_c;
}

void multigrid_restriction_gpu(double *d_val_c, double *d_val_f, subdomain *dm_c, subdomain *dm_f, int level)
{
    int nx_c, ny_c, nz_c;
    int i_gl_c, j_gl_c, k_gl_c;

    int nx_f = dm_f->nx;
    int ny_f = dm_f->ny;
    int nz_f = dm_f->nz;

    switch (aggregation_type)
    {
        case 0:
            nx_c = dm_c->nx;
            ny_c = dm_c->ny;
            nz_c = dm_c->nz;
            i_gl_c = 0;
            j_gl_c = 0;
            k_gl_c = 0;
            break;

        case 1:
            if (level == lv_aggregation - 1) {
                nx_c = dm_c->nx / comm_1d_x.nprocs;
                ny_c = dm_c->ny / comm_1d_y.nprocs;
                nz_c = dm_c->nz / comm_1d_z.nprocs;

                i_gl_c = nx_c * comm_1d_x.myrank;
                j_gl_c = ny_c * comm_1d_y.myrank;
                k_gl_c = nz_c * comm_1d_z.myrank;
            } else {
                nx_c = dm_c->nx;
                ny_c = dm_c->ny;
                nz_c = dm_c->nz;
                i_gl_c = j_gl_c = k_gl_c = 0;
            }
            break;

        case 2:
            nx_c = (level == lv_aggregation_x - 1)
                    ? dm_c->nx / comm_1d_x.nprocs
                    : dm_c->nx;
            i_gl_c = (level == lv_aggregation_x - 1)
                        ? nx_c * comm_1d_x.myrank
                        : 0;

            ny_c = (level == lv_aggregation_y - 1)
                    ? dm_c->ny / comm_1d_y.nprocs
                    : dm_c->ny;
            j_gl_c = (level == lv_aggregation_y - 1)
                        ? ny_c * comm_1d_y.myrank
                        : 0;

            nz_c = (level == lv_aggregation_z - 1)
                    ? dm_c->nz / comm_1d_z.nprocs
                    : dm_c->nz;
            k_gl_c = (level == lv_aggregation_z - 1)
                        ? nz_c * comm_1d_z.myrank
                        : 0;
            break;

        default:
            if (myrank == 0) {
                printf("[Error] Aggregation method should be 0, 1, or 2: %d\n",

                       aggregation_type);
            }
            MPI_Finalize();
            exit(EXIT_FAILURE);
    }

    multigrid_workspace_grow_gpu(&mg_restrict_dxf, &mg_restrict_dxf_capacity, (size_t)nx_f + 2);
    multigrid_workspace_grow_gpu(&mg_restrict_dyf, &mg_restrict_dyf_capacity, (size_t)ny_f + 2);
    multigrid_workspace_grow_gpu(&mg_restrict_dzf, &mg_restrict_dzf_capacity, (size_t)nz_f + 2);
    double *dxf = mg_restrict_dxf;
    double *dyf = mg_restrict_dyf;
    double *dzf = mg_restrict_dzf;

    int block1d = RESTRICTION_BLOCK_1D;
    build_dxf_kernel<<<(nx_f + 2 + block1d - 1) / block1d, block1d>>>(
        dxf, dm_c->dxm, dm_f->dxm, nx_f, level, lv_gdm_coarsest_x);
    build_dyf_kernel<<<(ny_f + 2 + block1d - 1) / block1d, block1d>>>(
        dyf, dm_c->dym, dm_f->dym, ny_f, level, lv_gdm_coarsest_y);
    build_dzf_kernel<<<(nz_f + 2 + block1d - 1) / block1d, block1d>>>(
        dzf, dm_c->dzm, dm_f->dzm, nz_f, level, lv_gdm_coarsest_z);

    CUDA_CHECK(cudaGetLastError());

    dim3 block(RESTRICTION_BLOCK_X, RESTRICTION_BLOCK_Y, RESTRICTION_BLOCK_Z);
    dim3 grid((nx_c + block.x - 1) / block.x,
              (ny_c + block.y - 1) / block.y,
              (nz_c + block.z - 1) / block.z);

    restriction_kernel<<<grid, block>>>(
        d_val_c, d_val_f,
        dxf, dyf, dzf,
        nx_c, ny_c, nz_c,
        nx_f, ny_f, nz_f,
        dm_c->ny, dm_c->nz,
        i_gl_c, j_gl_c, k_gl_c,
        level,
        lv_gdm_coarsest_x,
        lv_gdm_coarsest_y,
        lv_gdm_coarsest_z);

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

}

// __global__
// static void restriction_kernel(double *val_c, const double *val_f, const double *dxf, const double *dyf, const double *dzf,
//                                int nx_c, int ny_c, int nz_c,
//                                int nx_f, int ny_f, int nz_f,
//                                int dm_c_ny, int dm_c_nz,
//                                int level, int lv_gdm_coarsest_x, int lv_gdm_coarsest_y, int lv_gdm_coarsest_z)
// {
//     int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
//     int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
//     int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

//     if (i > nx_c || j > ny_c || k > nz_c) return;

//     int i_stride_f = (level >= lv_gdm_coarsest_x) ? 1 : 2;
//     int j_stride_f = (level >= lv_gdm_coarsest_y) ? 1 : 2;
//     int k_stride_f = (level >= lv_gdm_coarsest_z) ? 1 : 2;

//     int i_offset_f = (level >= lv_gdm_coarsest_x) ? 0 : 1;
//     int j_offset_f = (level >= lv_gdm_coarsest_y) ? 0 : 1;
//     int k_offset_f = (level >= lv_gdm_coarsest_z) ? 0 : 1;

//     // No MPI aggregation: local coarse index equals global coarse index.
//     int i_c = i;
//     int j_c = j;
//     int k_c = k;

//     int ip_f = i * i_stride_f;
//     int jp_f = j * j_stride_f;
//     int kp_f = k * k_stride_f;

//     int iz_f = ip_f - i_offset_f;
//     int jz_f = jp_f - j_offset_f;
//     int kz_f = kp_f - k_offset_f;

//     int ni_c = (dm_c_ny + 2) * (dm_c_nz + 2);
//     int nj_c = (dm_c_nz + 2);

//     int ni_f = (ny_f + 2) * (nz_f + 2);
//     int nj_f = (nz_f + 2);

//     double v000 = dxf[iz_f] * dyf[jz_f] * dzf[kz_f];
//     double v100 = dxf[ip_f] * dyf[jz_f] * dzf[kz_f];
//     double v010 = dxf[iz_f] * dyf[jp_f] * dzf[kz_f];
//     double v110 = dxf[ip_f] * dyf[jp_f] * dzf[kz_f];
//     double v001 = dxf[iz_f] * dyf[jz_f] * dzf[kp_f];
//     double v101 = dxf[ip_f] * dyf[jz_f] * dzf[kp_f];
//     double v011 = dxf[iz_f] * dyf[jp_f] * dzf[kp_f];
//     double v111 = dxf[ip_f] * dyf[jp_f] * dzf[kp_f];

//     double vol_c = v000 + v100 + v010 + v110 + v001 + v101 + v011 + v111;

//     val_c[IDX(i_c,j_c,k_c,ni_c,nj_c)] =
//         ( v000 * val_f[IDX(ip_f,jp_f,kp_f,ni_f,nj_f)]
//         + v100 * val_f[IDX(iz_f,jp_f,kp_f,ni_f,nj_f)]
//         + v010 * val_f[IDX(ip_f,jz_f,kp_f,ni_f,nj_f)]
//         + v110 * val_f[IDX(iz_f,jz_f,kp_f,ni_f,nj_f)]
//         + v001 * val_f[IDX(ip_f,jp_f,kz_f,ni_f,nj_f)]
//         + v101 * val_f[IDX(iz_f,jp_f,kz_f,ni_f,nj_f)]
//         + v011 * val_f[IDX(ip_f,jz_f,kz_f,ni_f,nj_f)]
//         + v111 * val_f[IDX(iz_f,jz_f,kz_f,ni_f,nj_f)] ) / vol_c;
// }

// void multigrid_restriction_gpu(double *d_val_c, double *d_val_f, subdomain *dm_c, subdomain *dm_f, int level)
// {
//     int nx_c = dm_c->nx;
//     int ny_c = dm_c->ny;
//     int nz_c = dm_c->nz;

//     int nx_f = dm_f->nx;
//     int ny_f = dm_f->ny;
//     int nz_f = dm_f->nz;

//     double *dxf = NULL;
//     double *dyf = NULL;
//     double *dzf = NULL;

//     CUDA_CHECK(cudaMalloc((void**)&dxf, (nx_f + 2) * sizeof(double)));
//     CUDA_CHECK(cudaMalloc((void**)&dyf, (ny_f + 2) * sizeof(double)));
//     CUDA_CHECK(cudaMalloc((void**)&dzf, (nz_f + 2) * sizeof(double)));

//     int block1d = RESTRICTION_BLOCK_1D;
//     build_dxf_kernel<<<(nx_f + 2 + block1d - 1) / block1d, block1d>>>(
//         dxf, dm_c->dxm, dm_f->dxm, nx_f, level, lv_gdm_coarsest_x);
//     build_dyf_kernel<<<(ny_f + 2 + block1d - 1) / block1d, block1d>>>(
//         dyf, dm_c->dym, dm_f->dym, ny_f, level, lv_gdm_coarsest_y);
//     build_dzf_kernel<<<(nz_f + 2 + block1d - 1) / block1d, block1d>>>(
//         dzf, dm_c->dzm, dm_f->dzm, nz_f, level, lv_gdm_coarsest_z);

//     CUDA_CHECK(cudaGetLastError());

//     dim3 block(RESTRICTION_BLOCK_X, RESTRICTION_BLOCK_Y, RESTRICTION_BLOCK_Z);
//     dim3 grid((nx_c + block.x - 1) / block.x,
//               (ny_c + block.y - 1) / block.y,
//               (nz_c + block.z - 1) / block.z);

//     restriction_kernel<<<grid, block>>>(
//         d_val_c, d_val_f,
//         dxf, dyf, dzf,
//         nx_c, ny_c, nz_c,
//         nx_f, ny_f, nz_f,
//         dm_c->ny, dm_c->nz,
//         level,
//         lv_gdm_coarsest_x,
//         lv_gdm_coarsest_y,
//         lv_gdm_coarsest_z);

//     CUDA_CHECK(cudaGetLastError());
//     CUDA_CHECK(cudaDeviceSynchronize());

//     CUDA_CHECK(cudaFree(dxf));
//     CUDA_CHECK(cudaFree(dyf));
//     CUDA_CHECK(cudaFree(dzf));
// }


#ifndef PROLONGATION_BLOCK_1D
#define PROLONGATION_BLOCK_1D 256
#endif

#ifndef PROLONGATION_BLOCK_X
#define PROLONGATION_BLOCK_X 8
#define PROLONGATION_BLOCK_Y 8
#define PROLONGATION_BLOCK_Z 4
#endif

__global__
static void build_dxp_dxn_kernel(double *dxp, double *dxn, const double *dxm_c, const double *dxm_f, int nx_c, int nx_f, int i_gl_c, int level, int lv_gdm_coarsest_x)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (level >= lv_gdm_coarsest_x) {
        // Same resolution in this direction.
        // In this case, fine index range is effectively the same as coarse.
        if (i <= nx_f + 1) {
            dxp[i] = dxm_c[i];
            dxn[i] = dxm_c[i];
        }
    } else {
        // One coarse cell corresponds to two fine cells.
        // One thread handles one coarse index i = 1..nx_c.
        if (i >= 1 && i <= nx_c) {
            int ip_f = 2 * i;

            dxp[ip_f - 1] = 0.5 * dxm_c[i - 1 + i_gl_c] + 0.5 * dxm_f[ip_f - 1];
            dxp[ip_f]     = 0.5 * dxm_c[i     + i_gl_c] - 0.5 * dxm_f[ip_f - 1];

            dxn[ip_f - 1] = 0.5 * dxm_c[i     + i_gl_c] - 0.5 * dxm_f[ip_f];
            dxn[ip_f]     = 0.5 * dxm_c[i + 1 + i_gl_c] + 0.5 * dxm_f[ip_f];
        }
    }
}

__global__
static void build_dyp_dyn_kernel(double *dyp, double *dyn, const double *dym_c, const double *dym_f, int ny_c, int ny_f, int j_gl_c, int level, int lv_gdm_coarsest_y)
{
    int j = blockIdx.x * blockDim.x + threadIdx.x;

    if (level >= lv_gdm_coarsest_y) {
        if (j <= ny_f + 1) {
            dyp[j] = dym_c[j];
            dyn[j] = dym_c[j];
        }
    } else {
        if (j >= 1 && j <= ny_c) {
            int jp_f = 2 * j;

            dyp[jp_f - 1] = 0.5 * dym_c[j - 1 + j_gl_c] + 0.5 * dym_f[jp_f - 1];
            dyp[jp_f]     = 0.5 * dym_c[j     + j_gl_c] - 0.5 * dym_f[jp_f - 1];

            dyn[jp_f - 1] = 0.5 * dym_c[j     + j_gl_c] - 0.5 * dym_f[jp_f];
            dyn[jp_f]     = 0.5 * dym_c[j + 1 + j_gl_c] + 0.5 * dym_f[jp_f];
        }
    }
}

__global__
static void build_dzp_dzn_kernel(double *dzp, double *dzn, const double *dzm_c, const double *dzm_f,  int nz_c, int nz_f, int k_gl_c, int level, int lv_gdm_coarsest_z)
{
    int k = blockIdx.x * blockDim.x + threadIdx.x;

    if (level >= lv_gdm_coarsest_z) {
        if (k <= nz_f + 1) {
            dzp[k] = dzm_c[k];
            dzn[k] = dzm_c[k];
        }
    } else {
        if (k >= 1 && k <= nz_c) {
            int kp_f = 2 * k;

            dzp[kp_f - 1] = 0.5 * dzm_c[k - 1 + k_gl_c] + 0.5 * dzm_f[kp_f - 1];
            dzp[kp_f]     = 0.5 * dzm_c[k     + k_gl_c] - 0.5 * dzm_f[kp_f - 1];

            dzn[kp_f - 1] = 0.5 * dzm_c[k     + k_gl_c] - 0.5 * dzm_f[kp_f];
            dzn[kp_f]     = 0.5 * dzm_c[k + 1 + k_gl_c] + 0.5 * dzm_f[kp_f];
        }
    }
}

__device__
static inline double prolong_interp8(const double *val_c,
                                     int ni_c, int nj_c,

                                     int iz_c, int jz_c, int kz_c,
                                     int im_c, int jm_c, int km_c,
                                     int ip_c, int jp_c, int kp_c,
                                     const double *wx0, const double *wx1,
                                     const double *wy0, const double *wy1,
                                     const double *wz0, const double *wz1,
                                     int ifine, int jfine, int kfine,
                                     int sx, int sy, int sz)
{
    // sx = 0: fine point is on the negative/left side in x, use dxp and neighbor im_c.
    // sx = 1: fine point is on the positive/right side in x, use dxn and neighbor ip_c.
    int x0 = iz_c;
    int x1 = (sx == 0) ? im_c : ip_c;
    double ax0 = (sx == 0) ? wx0[ifine] : wx0[ifine];
    double ax1 = (sx == 0) ? wx0[ifine + 1] : wx0[ifine - 1];

    int y0 = jz_c;
    int y1 = (sy == 0) ? jm_c : jp_c;
    double ay0 = (sy == 0) ? wy0[jfine] : wy0[jfine];
    double ay1 = (sy == 0) ? wy0[jfine + 1] : wy0[jfine - 1];

    int z0 = kz_c;
    int z1 = (sz == 0) ? km_c : kp_c;
    double az0 = (sz == 0) ? wz0[kfine] : wz0[kfine];
    double az1 = (sz == 0) ? wz0[kfine + 1] : wz0[kfine - 1];

    // The pointer arguments are chosen so that wx0/wy0/wz0 are already dxp/dyp/dzp or dxn/dyn/dzn.
    // For sx/sy/sz = 0, the neighboring weight is at ip_f/jp_f/kp_f.
    // For sx/sy/sz = 1, the neighboring weight is at iz_f/jz_f/kz_f.
    double denom = (ax0 + ax1) * (ay0 + ay1) * (az0 + az1);

    double v =
        ax0 * ay0 * az0 * val_c[IDX(x0,y0,z0,ni_c,nj_c)] +
        ax1 * ay0 * az0 * val_c[IDX(x1,y0,z0,ni_c,nj_c)] +
        ax0 * ay1 * az0 * val_c[IDX(x0,y1,z0,ni_c,nj_c)] +
        ax0 * ay0 * az1 * val_c[IDX(x0,y0,z1,ni_c,nj_c)] +
        ax1 * ay1 * az0 * val_c[IDX(x1,y1,z0,ni_c,nj_c)] +
        ax1 * ay0 * az1 * val_c[IDX(x1,y0,z1,ni_c,nj_c)] +
        ax0 * ay1 * az1 * val_c[IDX(x0,y1,z1,ni_c,nj_c)] +
        ax1 * ay1 * az1 * val_c[IDX(x1,y1,z1,ni_c,nj_c)];

    return v / denom;
}

__global__
static void prolongation_linear_nonuniform_kernel(double *val_f, const double *val_c,
                                                  const double *dxp, const double *dxn, const double *dyp, const double *dyn, const double *dzp, const double *dzn,
                                                  int nx_c, int ny_c, int nz_c,
                                                  int nx_f, int ny_f, int nz_f,
                                                  int dm_c_ny, int dm_c_nz,
                                                  int i_gl_c, int j_gl_c, int k_gl_c,
                                                  int level, int lv_gdm_coarsest_x, int lv_gdm_coarsest_y, int lv_gdm_coarsest_z)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

    if (i > nx_c || j > ny_c || k > nz_c) return;

    int i_stride_f = (level >= lv_gdm_coarsest_x) ? 1 : 2;
    int j_stride_f = (level >= lv_gdm_coarsest_y) ? 1 : 2;
    int k_stride_f = (level >= lv_gdm_coarsest_z) ? 1 : 2;

    int i_offset_f = (level >= lv_gdm_coarsest_x) ? 0 : 1;
    int j_offset_f = (level >= lv_gdm_coarsest_y) ? 0 : 1;
    int k_offset_f = (level >= lv_gdm_coarsest_z) ? 0 : 1;

    /* Coarse arrays use global aggregated coordinates at transition levels. */
    int iz_c = i + i_gl_c;
    int jz_c = j + j_gl_c;
    int kz_c = k + k_gl_c;

    int im_c = iz_c - i_offset_f;
    int jm_c = jz_c - j_offset_f;
    int km_c = kz_c - k_offset_f;

    int ip_c = iz_c + i_offset_f;
    int jp_c = jz_c + j_offset_f;
    int kp_c = kz_c + k_offset_f;

    int ip_f = i * i_stride_f;
    int jp_f = j * j_stride_f;
    int kp_f = k * k_stride_f;

    int iz_f = ip_f - i_offset_f;
    int jz_f = jp_f - j_offset_f;
    int kz_f = kp_f - k_offset_f;

    int ni_c = (dm_c_ny + 2) * (dm_c_nz + 2);
    int nj_c = (dm_c_nz + 2);

    int ni_f = (ny_f + 2) * (nz_f + 2);
    int nj_f = (nz_f + 2);

    // Same as the CPU code, but written compactly.
    // sx/sy/sz = 0 writes iz_f/jz_f/kz_f using dxp/dyp/dzp.
    // sx/sy/sz = 1 writes ip_f/jp_f/kp_f using dxn/dyn/dzn.
    val_f[IDX(iz_f,jz_f,kz_f,ni_f,nj_f)] =
        prolong_interp8(val_c, ni_c, nj_c,
                        iz_c,jz_c,kz_c, im_c,jm_c,km_c, ip_c,jp_c,kp_c,
                        dxp, dxp, dyp, dyp, dzp, dzp,
                        iz_f,jz_f,kz_f, 0,0,0);

    if (i_offset_f == 1) {
        val_f[IDX(ip_f,jz_f,kz_f,ni_f,nj_f)] =
            prolong_interp8(val_c, ni_c, nj_c,
                            iz_c,jz_c,kz_c, im_c,jm_c,km_c, ip_c,jp_c,kp_c,
                            dxn, dxn, dyp, dyp, dzp, dzp,
                            ip_f,jz_f,kz_f, 1,0,0);
    }

    if (j_offset_f == 1) {
        val_f[IDX(iz_f,jp_f,kz_f,ni_f,nj_f)] =
            prolong_interp8(val_c, ni_c, nj_c,
                            iz_c,jz_c,kz_c, im_c,jm_c,km_c, ip_c,jp_c,kp_c,
                            dxp, dxp, dyn, dyn, dzp, dzp,
                            iz_f,jp_f,kz_f, 0,1,0);
    }

    if (i_offset_f == 1 && j_offset_f == 1) {
        val_f[IDX(ip_f,jp_f,kz_f,ni_f,nj_f)] =
            prolong_interp8(val_c, ni_c, nj_c,
                            iz_c,jz_c,kz_c, im_c,jm_c,km_c, ip_c,jp_c,kp_c,
                            dxn, dxn, dyn, dyn, dzp, dzp,
                            ip_f,jp_f,kz_f, 1,1,0);
    }

    if (k_offset_f == 1) {
        val_f[IDX(iz_f,jz_f,kp_f,ni_f,nj_f)] =
            prolong_interp8(val_c, ni_c, nj_c,
                            iz_c,jz_c,kz_c, im_c,jm_c,km_c, ip_c,jp_c,kp_c,
                            dxp, dxp, dyp, dyp, dzn, dzn,
                            iz_f,jz_f,kp_f, 0,0,1);
    }

    if (i_offset_f == 1 && k_offset_f == 1) {
        val_f[IDX(ip_f,jz_f,kp_f,ni_f,nj_f)] =
            prolong_interp8(val_c, ni_c, nj_c,
                            iz_c,jz_c,kz_c, im_c,jm_c,km_c, ip_c,jp_c,kp_c,
                            dxn, dxn, dyp, dyp, dzn, dzn,
                            ip_f,jz_f,kp_f, 1,0,1);
    }

    if (j_offset_f == 1 && k_offset_f == 1) {
        val_f[IDX(iz_f,jp_f,kp_f,ni_f,nj_f)] =
            prolong_interp8(val_c, ni_c, nj_c,
                            iz_c,jz_c,kz_c, im_c,jm_c,km_c, ip_c,jp_c,kp_c,
                            dxp, dxp, dyn, dyn, dzn, dzn,
                            iz_f,jp_f,kp_f, 0,1,1);
    }

    if (i_offset_f == 1 && j_offset_f == 1 && k_offset_f == 1) {
        val_f[IDX(ip_f,jp_f,kp_f,ni_f,nj_f)] =
            prolong_interp8(val_c, ni_c, nj_c,
                            iz_c,jz_c,kz_c, im_c,jm_c,km_c, ip_c,jp_c,kp_c,
                            dxn, dxn, dyn, dyn, dzn, dzn,
                            ip_f,jp_f,kp_f, 1,1,1);
    }
}

void multigrid_prolongation_linear_on_nonuniform_grid_gpu(double *d_val_f, const double *d_val_c, subdomain *dm_f, subdomain *dm_c, int level)
{
    int nx_c, ny_c, nz_c;
    int i_gl_c, j_gl_c, k_gl_c;

    int nx_f = dm_f->nx;
    int ny_f = dm_f->ny;
    int nz_f = dm_f->nz;

    switch (aggregation_type)
    {
        case 0:
            nx_c = dm_c->nx;
            ny_c = dm_c->ny;
            nz_c = dm_c->nz;
            i_gl_c = 0;
            j_gl_c = 0;
            k_gl_c = 0;
            break;

        case 1:
            if (level == lv_aggregation - 1) {
                nx_c = dm_c->nx / comm_1d_x.nprocs;
                ny_c = dm_c->ny / comm_1d_y.nprocs;
                nz_c = dm_c->nz / comm_1d_z.nprocs;

                i_gl_c = nx_c * comm_1d_x.myrank;
                j_gl_c = ny_c * comm_1d_y.myrank;
                k_gl_c = nz_c * comm_1d_z.myrank;
            } else {
                nx_c = dm_c->nx;
                ny_c = dm_c->ny;
                nz_c = dm_c->nz;
                i_gl_c = j_gl_c = k_gl_c = 0;
            }
            break;

        case 2:
            nx_c = (level == lv_aggregation_x - 1)
                    ? dm_c->nx / comm_1d_x.nprocs
                    : dm_c->nx;
            i_gl_c = (level == lv_aggregation_x - 1)
                        ? nx_c * comm_1d_x.myrank
                        : 0;

            ny_c = (level == lv_aggregation_y - 1)
                    ? dm_c->ny / comm_1d_y.nprocs
                    : dm_c->ny;
            j_gl_c = (level == lv_aggregation_y - 1)
                        ? ny_c * comm_1d_y.myrank
                        : 0;

            nz_c = (level == lv_aggregation_z - 1)
                    ? dm_c->nz / comm_1d_z.nprocs
                    : dm_c->nz;
            k_gl_c = (level == lv_aggregation_z - 1)
                        ? nz_c * comm_1d_z.myrank
                        : 0;
            break;

        default:
            if (myrank == 0) {
                printf("[Error] Aggregation method should be 0, 1, or 2: %d\n",
                       aggregation_type);
            }
            MPI_Finalize();
            exit(EXIT_FAILURE);
    }


    const size_t x_required = (size_t)nx_f + 2;
    const size_t y_required = (size_t)ny_f + 2;
    const size_t z_required = (size_t)nz_f + 2;

    if (x_required > mg_prolong_x_capacity) {
        if (mg_prolong_dxp) CUDA_CHECK(cudaFree(mg_prolong_dxp));
        if (mg_prolong_dxn) CUDA_CHECK(cudaFree(mg_prolong_dxn));
        CUDA_CHECK(cudaMalloc((void **)&mg_prolong_dxp, x_required * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void **)&mg_prolong_dxn, x_required * sizeof(double)));
        mg_prolong_x_capacity = x_required;
    }
    if (y_required > mg_prolong_y_capacity) {
        if (mg_prolong_dyp) CUDA_CHECK(cudaFree(mg_prolong_dyp));
        if (mg_prolong_dyn) CUDA_CHECK(cudaFree(mg_prolong_dyn));
        CUDA_CHECK(cudaMalloc((void **)&mg_prolong_dyp, y_required * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void **)&mg_prolong_dyn, y_required * sizeof(double)));
        mg_prolong_y_capacity = y_required;
    }
    if (z_required > mg_prolong_z_capacity) {
        if (mg_prolong_dzp) CUDA_CHECK(cudaFree(mg_prolong_dzp));
        if (mg_prolong_dzn) CUDA_CHECK(cudaFree(mg_prolong_dzn));
        CUDA_CHECK(cudaMalloc((void **)&mg_prolong_dzp, z_required * sizeof(double)));
        CUDA_CHECK(cudaMalloc((void **)&mg_prolong_dzn, z_required * sizeof(double)));
        mg_prolong_z_capacity = z_required;
    }

    double *dxp = mg_prolong_dxp, *dxn = mg_prolong_dxn;
    double *dyp = mg_prolong_dyp, *dyn = mg_prolong_dyn;
    double *dzp = mg_prolong_dzp, *dzn = mg_prolong_dzn;

    int block1d = PROLONGATION_BLOCK_1D;

    build_dxp_dxn_kernel<<<(nx_f + 2 + block1d - 1) / block1d, block1d>>>(
        dxp, dxn, dm_c->dxm, dm_f->dxm,
        nx_c, nx_f, i_gl_c, level, lv_gdm_coarsest_x);

    build_dyp_dyn_kernel<<<(ny_f + 2 + block1d - 1) / block1d, block1d>>>(
        dyp, dyn, dm_c->dym, dm_f->dym,
        ny_c, ny_f, j_gl_c, level, lv_gdm_coarsest_y);


    build_dzp_dzn_kernel<<<(nz_f + 2 + block1d - 1) / block1d, block1d>>>(
        dzp, dzn, dm_c->dzm, dm_f->dzm,
        nz_c, nz_f, k_gl_c, level, lv_gdm_coarsest_z);

    CUDA_CHECK(cudaGetLastError());

    dim3 block(PROLONGATION_BLOCK_X, PROLONGATION_BLOCK_Y, PROLONGATION_BLOCK_Z);

    dim3 grid((nx_c + block.x - 1) / block.x,
              (ny_c + block.y - 1) / block.y,
              (nz_c + block.z - 1) / block.z);

    prolongation_linear_nonuniform_kernel<<<grid, block>>>( d_val_f, d_val_c,
        dxp, dxn, dyp, dyn, dzp, dzn,
        nx_c, ny_c, nz_c, nx_f, ny_f, nz_f,
        dm_c->ny, dm_c->nz,
        i_gl_c, j_gl_c, k_gl_c,
        level, lv_gdm_coarsest_x, lv_gdm_coarsest_y, lv_gdm_coarsest_z);

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

}

// __global__
// static void build_dxp_dxn_kernel(double *dxp, double *dxn, const double *dxm_c, const double *dxm_f, int nx_c, int nx_f, int level, int lv_gdm_coarsest_x)
// {
//     int i = blockIdx.x * blockDim.x + threadIdx.x;

//     if (level >= lv_gdm_coarsest_x) {
//         // Same resolution in this direction.
//         // In this case, fine index range is effectively the same as coarse.
//         if (i <= nx_f + 1) {
//             dxp[i] = dxm_c[i];
//             dxn[i] = dxm_c[i];
//         }
//     } else {
//         // One coarse cell corresponds to two fine cells.
//         // One thread handles one coarse index i = 1..nx_c.
//         if (i >= 1 && i <= nx_c) {
//             int ip_f = 2 * i;

//             dxp[ip_f - 1] = 0.5 * dxm_c[i - 1] + 0.5 * dxm_f[ip_f - 1];
//             dxp[ip_f]     = 0.5 * dxm_c[i]     - 0.5 * dxm_f[ip_f - 1];

//             dxn[ip_f - 1] = 0.5 * dxm_c[i]     - 0.5 * dxm_f[ip_f];
//             dxn[ip_f]     = 0.5 * dxm_c[i + 1] + 0.5 * dxm_f[ip_f];
//         }
//     }
// }

// __global__
// static void build_dyp_dyn_kernel(double *dyp, double *dyn, const double *dym_c, const double *dym_f, int ny_c, int ny_f, int level, int lv_gdm_coarsest_y)
// {
//     int j = blockIdx.x * blockDim.x + threadIdx.x;

//     if (level >= lv_gdm_coarsest_y) {
//         if (j <= ny_f + 1) {
//             dyp[j] = dym_c[j];
//             dyn[j] = dym_c[j];
//         }
//     } else {
//         if (j >= 1 && j <= ny_c) {
//             int jp_f = 2 * j;

//             dyp[jp_f - 1] = 0.5 * dym_c[j - 1] + 0.5 * dym_f[jp_f - 1];
//             dyp[jp_f]     = 0.5 * dym_c[j]     - 0.5 * dym_f[jp_f - 1];

//             dyn[jp_f - 1] = 0.5 * dym_c[j]     - 0.5 * dym_f[jp_f];
//             dyn[jp_f]     = 0.5 * dym_c[j + 1] + 0.5 * dym_f[jp_f];
//         }
//     }
// }

// __global__
// static void build_dzp_dzn_kernel(double *dzp, double *dzn, const double *dzm_c, const double *dzm_f,  int nz_c, int nz_f, int level, int lv_gdm_coarsest_z)
// {
//     int k = blockIdx.x * blockDim.x + threadIdx.x;

//     if (level >= lv_gdm_coarsest_z) {
//         if (k <= nz_f + 1) {
//             dzp[k] = dzm_c[k];
//             dzn[k] = dzm_c[k];
//         }
//     } else {
//         if (k >= 1 && k <= nz_c) {
//             int kp_f = 2 * k;

//             dzp[kp_f - 1] = 0.5 * dzm_c[k - 1] + 0.5 * dzm_f[kp_f - 1];
//             dzp[kp_f]     = 0.5 * dzm_c[k]     - 0.5 * dzm_f[kp_f - 1];

//             dzn[kp_f - 1] = 0.5 * dzm_c[k]     - 0.5 * dzm_f[kp_f];
//             dzn[kp_f]     = 0.5 * dzm_c[k + 1] + 0.5 * dzm_f[kp_f];
//         }
//     }
// }

// __device__
// static inline double prolong_interp8(const double *val_c,
//                                      int ni_c, int nj_c,
//                                      int iz_c, int jz_c, int kz_c,
//                                      int im_c, int jm_c, int km_c,
//                                      int ip_c, int jp_c, int kp_c,
//                                      const double *wx0, const double *wx1,
//                                      const double *wy0, const double *wy1,
//                                      const double *wz0, const double *wz1,
//                                      int ifine, int jfine, int kfine,
//                                      int sx, int sy, int sz)
// {
//     // sx = 0: fine point is on the negative/left side in x, use dxp and neighbor im_c.
//     // sx = 1: fine point is on the positive/right side in x, use dxn and neighbor ip_c.
//     int x0 = iz_c;
//     int x1 = (sx == 0) ? im_c : ip_c;
//     double ax0 = (sx == 0) ? wx0[ifine] : wx0[ifine];
//     double ax1 = (sx == 0) ? wx0[ifine + 1] : wx0[ifine - 1];

//     int y0 = jz_c;
//     int y1 = (sy == 0) ? jm_c : jp_c;
//     double ay0 = (sy == 0) ? wy0[jfine] : wy0[jfine];
//     double ay1 = (sy == 0) ? wy0[jfine + 1] : wy0[jfine - 1];

//     int z0 = kz_c;
//     int z1 = (sz == 0) ? km_c : kp_c;
//     double az0 = (sz == 0) ? wz0[kfine] : wz0[kfine];
//     double az1 = (sz == 0) ? wz0[kfine + 1] : wz0[kfine - 1];

//     // The pointer arguments are chosen so that wx0/wy0/wz0 are already dxp/dyp/dzp or dxn/dyn/dzn.
//     // For sx/sy/sz = 0, the neighboring weight is at ip_f/jp_f/kp_f.
//     // For sx/sy/sz = 1, the neighboring weight is at iz_f/jz_f/kz_f.
//     double denom = (ax0 + ax1) * (ay0 + ay1) * (az0 + az1);

//     double v =
//         ax0 * ay0 * az0 * val_c[IDX(x0,y0,z0,ni_c,nj_c)] +
//         ax1 * ay0 * az0 * val_c[IDX(x1,y0,z0,ni_c,nj_c)] +
//         ax0 * ay1 * az0 * val_c[IDX(x0,y1,z0,ni_c,nj_c)] +
//         ax0 * ay0 * az1 * val_c[IDX(x0,y0,z1,ni_c,nj_c)] +
//         ax1 * ay1 * az0 * val_c[IDX(x1,y1,z0,ni_c,nj_c)] +
//         ax1 * ay0 * az1 * val_c[IDX(x1,y0,z1,ni_c,nj_c)] +
//         ax0 * ay1 * az1 * val_c[IDX(x0,y1,z1,ni_c,nj_c)] +
//         ax1 * ay1 * az1 * val_c[IDX(x1,y1,z1,ni_c,nj_c)];

//     return v / denom;
// }

// __global__
// static void prolongation_linear_nonuniform_kernel(double *val_f, const double *val_c,
//                                                   const double *dxp, const double *dxn, const double *dyp, const double *dyn, const double *dzp, const double *dzn,
//                                                   int nx_c, int ny_c, int nz_c,
//                                                   int nx_f, int ny_f, int nz_f,
//                                                   int dm_c_ny, int dm_c_nz,
//                                                   int level, int lv_gdm_coarsest_x, int lv_gdm_coarsest_y, int lv_gdm_coarsest_z)
// {
//     int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
//     int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
//     int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

//     if (i > nx_c || j > ny_c || k > nz_c) return;

//     int i_stride_f = (level >= lv_gdm_coarsest_x) ? 1 : 2;
//     int j_stride_f = (level >= lv_gdm_coarsest_y) ? 1 : 2;
//     int k_stride_f = (level >= lv_gdm_coarsest_z) ? 1 : 2;

//     int i_offset_f = (level >= lv_gdm_coarsest_x) ? 0 : 1;
//     int j_offset_f = (level >= lv_gdm_coarsest_y) ? 0 : 1;
//     int k_offset_f = (level >= lv_gdm_coarsest_z) ? 0 : 1;

//     // No MPI aggregation.
//     int iz_c = i;
//     int jz_c = j;
//     int kz_c = k;

//     int im_c = iz_c - i_offset_f;
//     int jm_c = jz_c - j_offset_f;
//     int km_c = kz_c - k_offset_f;

//     int ip_c = iz_c + i_offset_f;
//     int jp_c = jz_c + j_offset_f;
//     int kp_c = kz_c + k_offset_f;

//     int ip_f = i * i_stride_f;
//     int jp_f = j * j_stride_f;
//     int kp_f = k * k_stride_f;

//     int iz_f = ip_f - i_offset_f;
//     int jz_f = jp_f - j_offset_f;
//     int kz_f = kp_f - k_offset_f;

//     int ni_c = (dm_c_ny + 2) * (dm_c_nz + 2);
//     int nj_c = (dm_c_nz + 2);

//     int ni_f = (ny_f + 2) * (nz_f + 2);
//     int nj_f = (nz_f + 2);

//     // Same as the CPU code, but written compactly.
//     // sx/sy/sz = 0 writes iz_f/jz_f/kz_f using dxp/dyp/dzp.
//     // sx/sy/sz = 1 writes ip_f/jp_f/kp_f using dxn/dyn/dzn.
//     val_f[IDX(iz_f,jz_f,kz_f,ni_f,nj_f)] =
//         prolong_interp8(val_c, ni_c, nj_c,
//                         iz_c,jz_c,kz_c, im_c,jm_c,km_c, ip_c,jp_c,kp_c,
//                         dxp, dxp, dyp, dyp, dzp, dzp,
//                         iz_f,jz_f,kz_f, 0,0,0);

//     if (i_offset_f == 1) {
//         val_f[IDX(ip_f,jz_f,kz_f,ni_f,nj_f)] =
//             prolong_interp8(val_c, ni_c, nj_c,
//                             iz_c,jz_c,kz_c, im_c,jm_c,km_c, ip_c,jp_c,kp_c,
//                             dxn, dxn, dyp, dyp, dzp, dzp,
//                             ip_f,jz_f,kz_f, 1,0,0);
//     }

//     if (j_offset_f == 1) {
//         val_f[IDX(iz_f,jp_f,kz_f,ni_f,nj_f)] =
//             prolong_interp8(val_c, ni_c, nj_c,
//                             iz_c,jz_c,kz_c, im_c,jm_c,km_c, ip_c,jp_c,kp_c,
//                             dxp, dxp, dyn, dyn, dzp, dzp,
//                             iz_f,jp_f,kz_f, 0,1,0);
//     }

//     if (i_offset_f == 1 && j_offset_f == 1) {
//         val_f[IDX(ip_f,jp_f,kz_f,ni_f,nj_f)] =
//             prolong_interp8(val_c, ni_c, nj_c,
//                             iz_c,jz_c,kz_c, im_c,jm_c,km_c, ip_c,jp_c,kp_c,
//                             dxn, dxn, dyn, dyn, dzp, dzp,
//                             ip_f,jp_f,kz_f, 1,1,0);
//     }

//     if (k_offset_f == 1) {
//         val_f[IDX(iz_f,jz_f,kp_f,ni_f,nj_f)] =
//             prolong_interp8(val_c, ni_c, nj_c,
//                             iz_c,jz_c,kz_c, im_c,jm_c,km_c, ip_c,jp_c,kp_c,
//                             dxp, dxp, dyp, dyp, dzn, dzn,
//                             iz_f,jz_f,kp_f, 0,0,1);
//     }

//     if (i_offset_f == 1 && k_offset_f == 1) {
//         val_f[IDX(ip_f,jz_f,kp_f,ni_f,nj_f)] =
//             prolong_interp8(val_c, ni_c, nj_c,
//                             iz_c,jz_c,kz_c, im_c,jm_c,km_c, ip_c,jp_c,kp_c,
//                             dxn, dxn, dyp, dyp, dzn, dzn,
//                             ip_f,jz_f,kp_f, 1,0,1);
//     }

//     if (j_offset_f == 1 && k_offset_f == 1) {

//         val_f[IDX(iz_f,jp_f,kp_f,ni_f,nj_f)] =
//             prolong_interp8(val_c, ni_c, nj_c,
//                             iz_c,jz_c,kz_c, im_c,jm_c,km_c, ip_c,jp_c,kp_c,
//                             dxp, dxp, dyn, dyn, dzn, dzn,
//                             iz_f,jp_f,kp_f, 0,1,1);
//     }

//     if (i_offset_f == 1 && j_offset_f == 1 && k_offset_f == 1) {
//         val_f[IDX(ip_f,jp_f,kp_f,ni_f,nj_f)] =
//             prolong_interp8(val_c, ni_c, nj_c,
//                             iz_c,jz_c,kz_c, im_c,jm_c,km_c, ip_c,jp_c,kp_c,
//                             dxn, dxn, dyn, dyn, dzn, dzn,
//                             ip_f,jp_f,kp_f, 1,1,1);
//     }
// }

// void multigrid_prolongation_linear_on_nonuniform_grid_gpu(double *d_val_f, const double *d_val_c, subdomain *dm_f, subdomain *dm_c, int level)
// {
//     int nx_c = dm_c->nx;
//     int ny_c = dm_c->ny;
//     int nz_c = dm_c->nz;

//     int nx_f = dm_f->nx;
//     int ny_f = dm_f->ny;
//     int nz_f = dm_f->nz;

//     double *dxp = NULL, *dxn = NULL;
//     double *dyp = NULL, *dyn = NULL;
//     double *dzp = NULL, *dzn = NULL;

//     CUDA_CHECK(cudaMalloc((void**)&dxp, (nx_f + 2) * sizeof(double)));
//     CUDA_CHECK(cudaMalloc((void**)&dxn, (nx_f + 2) * sizeof(double)));
//     CUDA_CHECK(cudaMalloc((void**)&dyp, (ny_f + 2) * sizeof(double)));
//     CUDA_CHECK(cudaMalloc((void**)&dyn, (ny_f + 2) * sizeof(double)));
//     CUDA_CHECK(cudaMalloc((void**)&dzp, (nz_f + 2) * sizeof(double)));
//     CUDA_CHECK(cudaMalloc((void**)&dzn, (nz_f + 2) * sizeof(double)));

//     int block1d = PROLONGATION_BLOCK_1D;

//     build_dxp_dxn_kernel<<<(nx_f + 2 + block1d - 1) / block1d, block1d>>>(
//         dxp, dxn, dm_c->dxm, dm_f->dxm,
//         nx_c, nx_f, level, lv_gdm_coarsest_x);

//     build_dyp_dyn_kernel<<<(ny_f + 2 + block1d - 1) / block1d, block1d>>>(
//         dyp, dyn, dm_c->dym, dm_f->dym,
//         ny_c, ny_f, level, lv_gdm_coarsest_y);

//     build_dzp_dzn_kernel<<<(nz_f + 2 + block1d - 1) / block1d, block1d>>>(
//         dzp, dzn, dm_c->dzm, dm_f->dzm,
//         nz_c, nz_f, level, lv_gdm_coarsest_z);

//     CUDA_CHECK(cudaGetLastError());

//     dim3 block(PROLONGATION_BLOCK_X, PROLONGATION_BLOCK_Y, PROLONGATION_BLOCK_Z);

//     dim3 grid((nx_c + block.x - 1) / block.x,
//               (ny_c + block.y - 1) / block.y,
//               (nz_c + block.z - 1) / block.z);

//     prolongation_linear_nonuniform_kernel<<<grid, block>>>( d_val_f, d_val_c,
//         dxp, dxn, dyp, dyn, dzp, dzn,
//         nx_c, ny_c, nz_c, nx_f, ny_f, nz_f,
//         dm_c->ny, dm_c->nz,
//         level, lv_gdm_coarsest_x, lv_gdm_coarsest_y, lv_gdm_coarsest_z);

//     CUDA_CHECK(cudaGetLastError());
//     CUDA_CHECK(cudaDeviceSynchronize());

//     CUDA_CHECK(cudaFree(dxp));
//     CUDA_CHECK(cudaFree(dxn));
//     CUDA_CHECK(cudaFree(dyp));
//     CUDA_CHECK(cudaFree(dyn));
//     CUDA_CHECK(cudaFree(dzp));
//     CUDA_CHECK(cudaFree(dzn));
// }


// 汝?ε댆瓦쇾틙




__global__
void residual_init_kernel(double *rsd, int nx, int ny, int nz)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    int j = blockIdx.y * blockDim.y + threadIdx.y;

    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if(i <= nx+1 && j <= ny+1 && k <= nz+1)
    {
        int ni=(ny+2)*(nz+2);
        int nj=(nz+2);

        int idx = i*ni + j*nj + k;

        rsd[idx]=0.0;
    }
}

__global__
void residual_update_kernel(double *rsd,  const double *rhs, int nx, int ny, int nz)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;

    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;

    int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

    if(i<=nx && j<=ny && k<=nz)
    {
        int ni=(ny+2)*(nz+2);
        int nj=(nz+2);

        int idx = i*ni + j*nj + k;

        rsd[idx] = rhs[idx] - rsd[idx];
    }
}

void multigrid_residual_gpu(double *d_rsd, double *d_coef, double *d_x, double *d_rhs, subdomain *dm, int is_aggregated[3])
{
    int nx = dm->nx;
    int ny = dm->ny;
    int nz = dm->nz;

    dim3 block(8,8,8);

    dim3 grid((nx+2 + block.x-1)/block.x, (ny+2 + block.y-1)/block.y, (nz+2 + block.z-1)/block.z);

    residual_init_kernel<<<grid,block>>>(d_rsd, nx, ny, nz);

    CUDA_CHECK(cudaGetLastError());

    /*
       d_rsd = A*x
    */
    mv_mul_poisson_matrix_gpu(d_rsd, d_coef, d_x, dm, is_aggregated);

    dim3 grid_inner((nx + block.x-1)/block.x, (ny + block.y-1)/block.y, (nz + block.z-1)/block.z);

    residual_update_kernel<<<grid_inner,block>>>(d_rsd, d_rhs, nx, ny, nz);

    CUDA_CHECK(cudaGetLastError());

    CUDA_CHECK(cudaDeviceSynchronize());
}


__global__
static void set_zero_kernel(double *x, size_t size)
{
    size_t idx = (size_t)blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < size) {
        x[idx] = 0.0;
    }
}

__global__
static void coarsest_direct_solve_kernel(double *x, const double *rhs, const double *coef, int nx, int ny, int nz)
{
    int ni = (ny + 2) * (nz + 2);
    int nj = (nz + 2);

    /*
       This kernel is only for nx*ny*nz == 1.
       The only interior point is (1,1,1).
    */
    if (blockIdx.x == 0 && threadIdx.x == 0) {
        int idx = IDX(1, 1, 1, ni, nj);

        int coef_size_one = (nx + 1) * (ny + 1) * (nz + 1);

        int coef_idx =
            0 * coef_size_one
          + 1 * (ny + 1) * (nz + 1)
          + 1 * (nz + 1)
          + 1;

        x[idx] = rhs[idx] / coef[coef_idx];
    }
}

void multigrid_solve_coarset_level_gpu(double *x, matrix_poisson *a_poisson, double *rhs, subdomain *dm, int maxiteration, double tolerance, double omega, int is_aggregated[3])
{
    int nx = dm->nx;
    int ny = dm->ny;
    int nz = dm->nz;

    size_t size = (size_t)(nx + 2) * (size_t)(ny + 2) * (size_t)(nz + 2);

    int block1d = 256;
    int grid1d = (int)((size + block1d - 1) / block1d);

    set_zero_kernel<<<grid1d, block1d>>>(x, size);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // if (myrank == 0) {
    //     printf("[MG] Obtain the solution on a grid in the coarsest level\n");
    // }

    if ((nx * ny * nz == 1) && is_aggregated[0] && is_aggregated[1] && is_aggregated[2])
    {
        coarsest_direct_solve_kernel<<<1,1>>>(x, rhs, a_poisson->coeff, nx, ny, nz);

        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
    }
    else
    {
        rbgs_solver_poisson_matrix_gpu(x, a_poisson->coeff, rhs, dm, maxiteration, tolerance, omega, is_aggregated);
    }



    
}

























void multigrid_solve_vcycle_gpu(double *d_sol, matrix_poisson *a_poisson, double *d_rhs, subdomain *sdm, int maxiteration, double tolerance, double omega_sor)
{
    int rank = -1;
    int mpi_inited = 0;

    MPI_Initialized(&mpi_inited);
    if (mpi_inited) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    } else {
        rank = myrank;
    }
    
    switch (aggregation_type)
    {
        case 0:
            multigrid_no_aggregation_vcycle_solver_gpu(d_sol, a_poisson, d_rhs, sdm, maxiteration, tolerance, omega_sor);
            break;

        case 1:
            multigrid_single_aggregation_vcycle_solver_gpu(d_sol, a_poisson, d_rhs, sdm, maxiteration, tolerance, omega_sor);
            break;

        case 2:
            multigrid_adaptive_aggregation_vcycle_solver_gpu(d_sol, a_poisson, d_rhs, sdm, maxiteration, tolerance, omega_sor);
            break;

        default:
            if (myrank == 0) {
                fprintf(stderr, "[Error] Aggregation method should be 0, 1, or 2. Got %d\n", aggregation_type);
            }
            exit(EXIT_FAILURE);
    }
}

__global__
static void vector_add_kernel(double *x, const double *corr, size_t n)
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        x[idx] += corr[idx];
    }
}

void multigrid_no_aggregation_vcycle_solver_gpu(double *sol, matrix_poisson *a_poisson, double *rhs, subdomain *sdm, int maxiteration, double tolerance, double omega_sor)
{
    int l, cyc;
    int nx = sdm->nx;
    int ny = sdm->ny;
    int nz = sdm->nz;

    // double rsd_val, res0tol;
    double rsd_local, rsd_val;
    double res0_local, res0tol;

    size_t size = (size_t)(nx + 2) * (size_t)(ny + 2) * (size_t)(nz + 2);

    double *rsd = NULL;
    CUDA_CHECK(cudaMalloc((void**)&rsd, size * sizeof(double)));
    CUDA_CHECK(cudaMemset(rsd, 0, size * sizeof(double)));

    // timer_stamp0(STAMP_COMP);
    // timer_stamp0(STAMP_COMM_NEIGHBOR);
    // timer_stamp0(STAMP_COMM_ALLREDUCE);
    // timer_stamp0(STAMP_residual);
    multigrid_residual_gpu(rsd, a_poisson->coeff, sol, rhs, sdm, sdm->is_aggregated);
    // timer_stamp(18,STAMP_residual);
    vv_dot_3d_matrix_gpu(&res0tol, rsd, rsd, nx, ny, nz, sdm->is_aggregated);

    for (cyc = 1; cyc <= n_vcycles; cyc++)
    {
        // timer_stamp0(STAMP_LEVEL);
        // timer_stamp0(STAMP_smooth);
        rbgs_iterator_poisson_matrix_gpu(sol, a_poisson->coeff, rhs, sdm, maxiteration, omega_sor, sdm->is_aggregated);
        // timer_stamp(15,STAMP_smooth);
        
        // timer_stamp0(STAMP_residual);
        multigrid_residual_gpu(rsd, a_poisson->coeff, sol, rhs, sdm, sdm->is_aggregated);
        //  timer_stamp(18,STAMP_residual);
        // timer_stamp(5,STAMP_LEVEL);

        // timer_stamp0(STAMP_restriction);
        multigrid_restriction_gpu(mg_sdm[1].b, rsd, &mg_sdm[1], sdm, 0);
        // timer_stamp(16,STAMP_restriction);

        // if (myrank == 0)
        //     printf("[MG] Restriction from level 0 to level 1\n");

        for (l = 1; l <= n_levels - 1; l++)
        {
            size_t level_size =
                (size_t)(mg_sdm[l].nx + 2) *
                (size_t)(mg_sdm[l].ny + 2) *
                (size_t)(mg_sdm[l].nz + 2);

            CUDA_CHECK(cudaMemset(mg_sdm[l].x, 0, level_size * sizeof(double)));

            // timer_stamp0(STAMP_LEVEL);
            // timer_stamp0(STAMP_smooth);
            rbgs_iterator_poisson_matrix_gpu( mg_sdm[l].x, mg_a_poisson[l].coeff, mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
            // timer_stamp(15,STAMP_smooth);

            // timer_stamp0(STAMP_residual);
            multigrid_residual_gpu(mg_sdm[l].r, mg_a_poisson[l].coeff, mg_sdm[l].x, mg_sdm[l].b, &mg_sdm[l], mg_sdm[l].is_aggregated);
            // timer_stamp(18,STAMP_residual);
            // timer_stamp(6,STAMP_LEVEL);

            // timer_stamp0(STAMP_restriction);
            multigrid_restriction_gpu(mg_sdm[l+1].b, mg_sdm[l].r, &mg_sdm[l+1], &mg_sdm[l],l);
            // timer_stamp(16,STAMP_restriction);

            // if (myrank == 0)
            //     printf("[MG] Restriction from level %d to level %d\n", l, l+1);
        }

        // timer_stamp0(STAMP_LEVEL);
        // timer_stamp0(STAMP_smooth);
        multigrid_solve_coarset_level_gpu(mg_sdm[n_levels].x, &mg_a_poisson[n_levels], mg_sdm[n_levels].b, &mg_sdm[n_levels], 1000, tolerance, omega_sor, mg_sdm[n_levels].is_aggregated);
        // timer_stamp(15,STAMP_smooth);

        // timer_stamp0(STAMP_residual);
        multigrid_residual_gpu(mg_sdm[n_levels].r, mg_a_poisson[n_levels].coeff, mg_sdm[n_levels].x, mg_sdm[n_levels].b, &mg_sdm[n_levels], mg_sdm[n_levels].is_aggregated);
        // timer_stamp(18,STAMP_residual);
        // timer_stamp(7,STAMP_LEVEL);

        if (myrank == 0) {
            int ni = (mg_sdm[n_levels].ny + 2) *
                     (mg_sdm[n_levels].nz + 2);
            int nj = mg_sdm[n_levels].nz + 2;
            int idx = IDX(1,1,1,ni,nj);

            double h_x111, h_r111;
            CUDA_CHECK(cudaMemcpy(&h_x111, &mg_sdm[n_levels].x[idx], sizeof(double), cudaMemcpyDeviceToHost));

            CUDA_CHECK(cudaMemcpy(&h_r111, &mg_sdm[n_levels].r[idx], sizeof(double), cudaMemcpyDeviceToHost));

            // printf("[MG] Solution in the coarsest level : x(1,1,1) = %18.10e, residue = %18.10e\n", h_x111, h_r111);
        }

        for (l = n_levels - 1; l >= 1; l--)
        {
            // timer_stamp0(STAMP_COMM_NEIGHBOR);
            geometry_halocell_update_selectively_gpu(mg_sdm[l+1].x, &mg_sdm[l+1], mg_sdm[l+1].is_aggregated);
            // timer_stamp(11,STAMP_COMM_NEIGHBOR);

            // timer_stamp0(STAMP_COMP);
            // timer_stamp0(STAMP_prolongation);
            multigrid_prolongation_linear_on_nonuniform_grid_gpu( mg_sdm[l].r, mg_sdm[l+1].x, &mg_sdm[l], &mg_sdm[l+1], l);
            // timer_stamp(17,STAMP_prolongation);
            // timer_stamp0(STAMP_LEVEL);

            size_t level_size =
                (size_t)(mg_sdm[l].nx + 2) *
                (size_t)(mg_sdm[l].ny + 2) *
                (size_t)(mg_sdm[l].nz + 2);

            int block1d = 256;
            int grid1d = (level_size + block1d - 1) / block1d;

            vector_add_kernel<<<grid1d, block1d>>>(mg_sdm[l].x, mg_sdm[l].r, level_size);

            CUDA_CHECK(cudaGetLastError());

            // timer_stamp(10,STAMP_COMP);
            // timer_stamp0(STAMP_smooth);
            rbgs_iterator_poisson_matrix_gpu(mg_sdm[l].x, mg_a_poisson[l].coeff, mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
            // timer_stamp(15,STAMP_smooth);
            // timer_stamp(6,STAMP_LEVEL);
        }

        // timer_stamp0(STAMP_COMM_NEIGHBOR);
        geometry_halocell_update_selectively_gpu(mg_sdm[l+1].x, &mg_sdm[l+1], mg_sdm[l+1].is_aggregated);
        // timer_stamp(11,STAMP_COMM_NEIGHBOR);
        
        // timer_stamp0(STAMP_COMP);
        // timer_stamp0(STAMP_prolongation);
        multigrid_prolongation_linear_on_nonuniform_grid_gpu(rsd, mg_sdm[1].x, sdm, &mg_sdm[1], 0);
        // timer_stamp(17,STAMP_prolongation);

        // timer_stamp0(STAMP_LEVEL);

        int block1d = 256;
        int grid1d = (size + block1d - 1) / block1d;

        vector_add_kernel<<<grid1d, block1d>>>(sol, rsd, size);
        // timer_stamp(10,STAMP_COMP);

        CUDA_CHECK(cudaGetLastError());

        // timer_stamp0(STAMP_smooth);
        rbgs_iterator_poisson_matrix_gpu(sol, a_poisson->coeff, rhs, sdm, maxiteration, omega_sor, sdm->is_aggregated);
        // timer_stamp(15,STAMP_smooth);
        // timer_stamp(5,STAMP_LEVEL);

        // timer_stamp0(STAMP_residual);
        multigrid_residual_gpu(rsd, a_poisson->coeff, sol, rhs, sdm, sdm->is_aggregated);
        // timer_stamp(18,STAMP_residual);
        
        vv_dot_3d_matrix_gpu(&rsd_val, rsd, rsd, nx, ny, nz, sdm->is_aggregated);

        // if (myrank == 0) {
        //     printf("[MG] cycle = %d, Error = %e, Initial = %e, Relative = %e\n", cyc, sqrt(rsd_val), sqrt(res0tol), sqrt(rsd_val / res0tol));
        // }

        if (sqrt(rsd_val / res0tol) < tolerance)
            break;
    }

    // if (myrank == 0)
    //     printf("[MG] Total %d V-cycles end\n", cyc);

    CUDA_CHECK(cudaFree(rsd));
}

void multigrid_single_aggregation_vcycle_solver_gpu(double *sol, matrix_poisson *a_poisson, double *rhs, subdomain *sdm, int maxiteration, double tolerance, double omega_sor)
{
    int i, l, cyc;
    double rsd_val, res0tol;
    MPI_Datatype ddtype_temp, ddtype_send, ddtype_gatherv;
    MPI_Datatype ddtype_scatterv, ddtype_recv;
    int sizes[3], subsizes[3], starts[3];
    int r8size;
    MPI_Aint extent;

    if (lv_aggregation < 1 || lv_aggregation > n_levels) {
        if (myrank == 0)
            fprintf(stderr,
                    "[Error] Invalid single aggregation level %d (levels=%d)\n",
                    lv_aggregation, n_levels);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    int (*cart_coord)[3] =
        (int (*)[3])malloc((size_t)nprocs * sizeof(*cart_coord));
    int *cnt_gatherv = (int *)malloc((size_t)nprocs * sizeof(int));
    int *disps_gatherv = (int *)malloc((size_t)nprocs * sizeof(int));
    int *cnt_scatterv = (int *)malloc((size_t)nprocs * sizeof(int));
    int *disps_scatterv = (int *)malloc((size_t)nprocs * sizeof(int));
    if (!cart_coord || !cnt_gatherv || !disps_gatherv ||
        !cnt_scatterv || !disps_scatterv) {
        fprintf(stderr, "[rank %d] Failed to allocate aggregation metadata\n",
                myrank);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    for (i = 0; i < nprocs; ++i) {
        int ierr = MPI_Cart_coords(mpi_world_cart, i, 3, cart_coord[i]);
        if (ierr != MPI_SUCCESS) {
            fprintf(stderr, "[rank %d] MPI_Cart_coords failed\n", myrank);

            MPI_Abort(MPI_COMM_WORLD, ierr);
        }
    }

    const int nx_aggr = mg_sdm[lv_aggregation].nx;
    const int ny_aggr = mg_sdm[lv_aggregation].ny;
    const int nz_aggr = mg_sdm[lv_aggregation].nz;
    const int nx_local = nx_aggr / comm_1d_x.nprocs;
    const int ny_local = ny_aggr / comm_1d_y.nprocs;
    const int nz_local = nz_aggr / comm_1d_z.nprocs;
    const int ix = nx_local * comm_1d_x.myrank;
    const int iy = ny_local * comm_1d_y.myrank;
    const int iz = nz_local * comm_1d_z.myrank;

    sizes[0] = nx_aggr + 2;
    sizes[1] = ny_aggr + 2;
    sizes[2] = nz_aggr + 2;
    subsizes[0] = nx_local;
    subsizes[1] = ny_local;
    subsizes[2] = nz_local;
    starts[0] = ix + 1;
    starts[1] = iy + 1;
    starts[2] = iz + 1;
    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_temp);
    MPI_Type_size(MPI_DOUBLE, &r8size);
    extent = (MPI_Aint)r8size;
    MPI_Type_create_resized(ddtype_temp, 0, extent, &ddtype_send);
    MPI_Type_commit(&ddtype_send);
    MPI_Type_free(&ddtype_temp);

    starts[0] = starts[1] = starts[2] = 1;
    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_temp);
    MPI_Type_create_resized(ddtype_temp, 0, extent, &ddtype_gatherv);
    MPI_Type_commit(&ddtype_gatherv);
    MPI_Type_free(&ddtype_temp);

    for (i = 0; i < nprocs; ++i) {
        cnt_gatherv[i] = 1;
        disps_gatherv[i] =
            nz_local * cart_coord[i][2] +
            ny_local * cart_coord[i][1] * (nz_aggr + 2) +
            nx_local * cart_coord[i][0] * (ny_aggr + 2) * (nz_aggr + 2);
    }

    subsizes[0] = nx_local + 2;
    subsizes[1] = ny_local + 2;
    subsizes[2] = nz_local + 2;
    starts[0] = starts[1] = starts[2] = 0;
    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_temp);
    MPI_Type_create_resized(ddtype_temp, 0, extent, &ddtype_scatterv);
    MPI_Type_commit(&ddtype_scatterv);
    MPI_Type_free(&ddtype_temp);

    starts[0] = ix;
    starts[1] = iy;
    starts[2] = iz;
    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_temp);
    MPI_Type_create_resized(ddtype_temp, 0, extent, &ddtype_recv);
    MPI_Type_commit(&ddtype_recv);
    MPI_Type_free(&ddtype_temp);

    for (i = 0; i < nprocs; ++i) {
        cnt_scatterv[i] = 1;
        disps_scatterv[i] =
            nz_local * cart_coord[i][2] +
            ny_local * cart_coord[i][1] * (nz_aggr + 2) +
            nx_local * cart_coord[i][0] * (ny_aggr + 2) * (nz_aggr + 2);
    }

    const int nx = sdm->nx;
    const int ny = sdm->ny;
    const int nz = sdm->nz;
    const size_t fine_size =
        (size_t)(nx + 2) * (size_t)(ny + 2) * (size_t)(nz + 2);
    double *rsd = NULL;
    CUDA_CHECK(cudaMalloc((void **)&rsd, fine_size * sizeof(double)));
    CUDA_CHECK(cudaMemset(rsd, 0, fine_size * sizeof(double)));

    multigrid_residual_gpu(rsd, a_poisson->coeff, sol, rhs, sdm,
                           sdm->is_aggregated);
    vv_dot_3d_matrix_gpu(&res0tol, rsd, rsd, nx, ny, nz,
                         sdm->is_aggregated);

    for (cyc = 1; cyc <= n_vcycles; ++cyc) {
        rbgs_iterator_poisson_matrix_gpu(
            sol, a_poisson->coeff, rhs, sdm, maxiteration, omega_sor,
            sdm->is_aggregated);
        multigrid_residual_gpu(rsd, a_poisson->coeff, sol, rhs, sdm,
                               sdm->is_aggregated);
        multigrid_restriction_gpu(mg_sdm[1].b, rsd, &mg_sdm[1], sdm, 0);
        // if (myrank == 0)
        //     printf("[MG-GPU] Restriction from level 0 to level 1\n");

        /* All ranks descend independently until the aggregation level. */
        for (l = 1; l <= lv_aggregation - 1; ++l) {
            const size_t level_size =
                (size_t)(mg_sdm[l].nx + 2) *
                (size_t)(mg_sdm[l].ny + 2) *
                (size_t)(mg_sdm[l].nz + 2);
            CUDA_CHECK(cudaMemset(mg_sdm[l].x, 0,
                                  level_size * sizeof(double)));
            rbgs_iterator_poisson_matrix_gpu(
                mg_sdm[l].x, mg_a_poisson[l].coeff, mg_sdm[l].b,
                &mg_sdm[l], maxiteration, omega_sor,
                mg_sdm[l].is_aggregated);
            multigrid_residual_gpu(
                mg_sdm[l].r, mg_a_poisson[l].coeff, mg_sdm[l].x,
                mg_sdm[l].b, &mg_sdm[l], mg_sdm[l].is_aggregated);
            multigrid_restriction_gpu(
                mg_sdm[l + 1].b, mg_sdm[l].r, &mg_sdm[l + 1],
                &mg_sdm[l], l);
            // if (myrank == 0)
            //     printf("[MG-GPU] Restriction from level %d to level %d\n",
            //            l, l + 1);
        }

        CUDA_CHECK(cudaDeviceSynchronize());
        if (myrank == 0) {
            MPI_Gatherv(MPI_IN_PLACE, 0, ddtype_send,
                        mg_sdm[lv_aggregation].b, cnt_gatherv,
                        disps_gatherv, ddtype_gatherv, 0, MPI_COMM_WORLD);
        } else {
            MPI_Gatherv(mg_sdm[lv_aggregation].b, 1, ddtype_send,
                        mg_sdm[lv_aggregation].b, cnt_gatherv,
                        disps_gatherv, ddtype_gatherv, 0, MPI_COMM_WORLD);
        }
        CUDA_CHECK(cudaDeviceSynchronize());

        /* The aggregated hierarchy exists on every rank, but rank 0 owns solve. */
        if (myrank == 0) {
            for (l = lv_aggregation; l < n_levels; ++l) {
                const size_t level_size =
                    (size_t)(mg_sdm[l].nx + 2) *
                    (size_t)(mg_sdm[l].ny + 2) *
                    (size_t)(mg_sdm[l].nz + 2);
                CUDA_CHECK(cudaMemset(mg_sdm[l].x, 0,
                                      level_size * sizeof(double)));
                rbgs_iterator_poisson_matrix_gpu(
                    mg_sdm[l].x, mg_a_poisson[l].coeff, mg_sdm[l].b,
                    &mg_sdm[l], maxiteration, omega_sor,
                    mg_sdm[l].is_aggregated);
                multigrid_residual_gpu(
                    mg_sdm[l].r, mg_a_poisson[l].coeff, mg_sdm[l].x,
                    mg_sdm[l].b, &mg_sdm[l], mg_sdm[l].is_aggregated);
                multigrid_restriction_gpu(
                    mg_sdm[l + 1].b, mg_sdm[l].r, &mg_sdm[l + 1],
                    &mg_sdm[l], l);
                // printf("[MG-GPU] Restriction from level %d to level %d\n",
                //        l, l + 1);
            }

            multigrid_solve_coarset_level_gpu(
                mg_sdm[n_levels].x, &mg_a_poisson[n_levels],
                mg_sdm[n_levels].b, &mg_sdm[n_levels], 1000, tolerance,
                omega_sor, mg_sdm[n_levels].is_aggregated);
            multigrid_residual_gpu(
                mg_sdm[n_levels].r, mg_a_poisson[n_levels].coeff,
                mg_sdm[n_levels].x, mg_sdm[n_levels].b,
                &mg_sdm[n_levels], mg_sdm[n_levels].is_aggregated);

            const int ni = (mg_sdm[n_levels].ny + 2) *
                           (mg_sdm[n_levels].nz + 2);
            const int nj = mg_sdm[n_levels].nz + 2;
            const int idx = IDX(1, 1, 1, ni, nj);
            double h_x111, h_r111;
            CUDA_CHECK(cudaMemcpy(&h_x111, mg_sdm[n_levels].x + idx,
                                  sizeof(double), cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(&h_r111, mg_sdm[n_levels].r + idx,
                                  sizeof(double), cudaMemcpyDeviceToHost));
            // printf("[MG-GPU] Coarsest x(1,1,1)=%18.10e, residual=%18.10e\n",
            //        h_x111, h_r111);

            for (l = n_levels - 1; l >= lv_aggregation; --l) {
                multigrid_prolongation_linear_on_nonuniform_grid_gpu(
                    mg_sdm[l].r, mg_sdm[l + 1].x, &mg_sdm[l],
                    &mg_sdm[l + 1], l);
                const size_t level_size =
                    (size_t)(mg_sdm[l].nx + 2) *
                    (size_t)(mg_sdm[l].ny + 2) *
                    (size_t)(mg_sdm[l].nz + 2);
                const int block1d = 256;
                const int grid1d =
                    (int)((level_size + block1d - 1) / block1d);
                vector_add_kernel<<<grid1d, block1d>>>(
                    mg_sdm[l].x, mg_sdm[l].r, level_size);
                CUDA_CHECK(cudaGetLastError());
                rbgs_iterator_poisson_matrix_gpu(
                    mg_sdm[l].x, mg_a_poisson[l].coeff, mg_sdm[l].b,
                    &mg_sdm[l], maxiteration, omega_sor,
                    mg_sdm[l].is_aggregated);
            }
        }

        CUDA_CHECK(cudaDeviceSynchronize());
        if (myrank == 0) {
            MPI_Scatterv(mg_sdm[lv_aggregation].x, cnt_scatterv,
                         disps_scatterv, ddtype_scatterv, MPI_IN_PLACE, 0,
                         ddtype_recv, 0, MPI_COMM_WORLD);
        } else {
            MPI_Scatterv(mg_sdm[lv_aggregation].x, cnt_scatterv,
                         disps_scatterv, ddtype_scatterv,
                         mg_sdm[lv_aggregation].x, 1, ddtype_recv, 0,
                         MPI_COMM_WORLD);
        }
        CUDA_CHECK(cudaDeviceSynchronize());

        for (l = lv_aggregation - 1; l >= 1; --l) {
            geometry_halocell_update_selectively_gpu(
                mg_sdm[l + 1].x, &mg_sdm[l + 1],
                mg_sdm[l + 1].is_aggregated);
            multigrid_prolongation_linear_on_nonuniform_grid_gpu(
                mg_sdm[l].r, mg_sdm[l + 1].x, &mg_sdm[l],
                &mg_sdm[l + 1], l);
            const size_t level_size =
                (size_t)(mg_sdm[l].nx + 2) *
                (size_t)(mg_sdm[l].ny + 2) *
                (size_t)(mg_sdm[l].nz + 2);
            const int block1d = 256;
            const int grid1d =
                (int)((level_size + block1d - 1) / block1d);
            vector_add_kernel<<<grid1d, block1d>>>(
                mg_sdm[l].x, mg_sdm[l].r, level_size);
            CUDA_CHECK(cudaGetLastError());
            rbgs_iterator_poisson_matrix_gpu(
                mg_sdm[l].x, mg_a_poisson[l].coeff, mg_sdm[l].b,
                &mg_sdm[l], maxiteration, omega_sor,
                mg_sdm[l].is_aggregated);
        }

        geometry_halocell_update_selectively_gpu(
            mg_sdm[1].x, &mg_sdm[1], mg_sdm[1].is_aggregated);
        multigrid_prolongation_linear_on_nonuniform_grid_gpu(
            rsd, mg_sdm[1].x, sdm, &mg_sdm[1], 0);
        const int block1d = 256;
        const int grid1d =
            (int)((fine_size + block1d - 1) / block1d);
        vector_add_kernel<<<grid1d, block1d>>>(sol, rsd, fine_size);
        CUDA_CHECK(cudaGetLastError());

        rbgs_iterator_poisson_matrix_gpu(
            sol, a_poisson->coeff, rhs, sdm, maxiteration, omega_sor,
            sdm->is_aggregated);
        multigrid_residual_gpu(rsd, a_poisson->coeff, sol, rhs, sdm,
                               sdm->is_aggregated);
        vv_dot_3d_matrix_gpu(&rsd_val, rsd, rsd, nx, ny, nz,
                             sdm->is_aggregated);


        // if (myrank == 0)
        //     printf("[MG-GPU] cycle=%d Error=%e Initial=%e Relative=%e\n",
        //            cyc, sqrt(rsd_val), sqrt(res0tol),
        //            sqrt(rsd_val / res0tol));
        if (sqrt(rsd_val / res0tol) < tolerance)
            break;
    }

    // if (myrank == 0)
    //     printf("[MG-GPU] Total %d V-cycles end\n", cyc);

    CUDA_CHECK(cudaFree(rsd));
    MPI_Type_free(&ddtype_recv);
    MPI_Type_free(&ddtype_scatterv);
    MPI_Type_free(&ddtype_gatherv);
    MPI_Type_free(&ddtype_send);
    free(disps_scatterv);
    free(cnt_scatterv);
    free(disps_gatherv);
    free(cnt_gatherv);
    free(cart_coord);
}

static void adaptive_vector_add_gpu(double *x, const double *corr, size_t n)
{
    const int block = 256;
    const int grid = (int)((n + block - 1) / block);
    vector_add_kernel<<<grid, block>>>(x, corr, n);
    CUDA_CHECK(cudaGetLastError());
}

static void adaptive_print_coarsest_gpu(int cyc)
{
    const int ni = (mg_sdm[n_levels].ny + 2) * (mg_sdm[n_levels].nz + 2);
    const int nj = mg_sdm[n_levels].nz + 2;
    const int idx = IDX(1, 1, 1, ni, nj);
    double hx, hr;
    CUDA_CHECK(cudaMemcpy(&hx, mg_sdm[n_levels].x + idx, sizeof(double),
                          cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(&hr, mg_sdm[n_levels].r + idx, sizeof(double),
                          cudaMemcpyDeviceToHost));
    printf("[MG-GPU] cycle=%d coarsest x(1,1,1)=%.10e residual=%.10e\n", cyc, hx, hr);
}












































void multigrid_adaptive_aggregation_vcycle_solver_gpu(double *sol, matrix_poisson *a_poisson, double *rhs, subdomain *sdm, int maxiteration, double tolerance, double omega_sor)
{
    int l, i, cyc,ni,nj;
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
    cart_coord = (int (*)[3])malloc((size_t)nprocs * sizeof(*cart_coord));
    int ierr;

    /* No direction needs aggregation (normally the single-process case). */
    if (lv_aggregation_x == 0 && lv_aggregation_y == 0 &&
        lv_aggregation_z == 0) {
        free(cart_coord);
        multigrid_no_aggregation_vcycle_solver_gpu(
            sol, a_poisson, rhs, sdm, maxiteration, tolerance, omega_sor);
        return;
    }

    
    for (i = 0; i < nprocs; i++) {
        ierr = MPI_Cart_coords(mpi_world_cart, i, 3, cart_coord[i]);
        if (ierr != MPI_SUCCESS) {
            fprintf(stderr, "MPI_Cart_coords error\n");
            MPI_Abort(MPI_COMM_WORLD, ierr);
        }
    }
    // cnt_gatherv_x = malloc((comm_1d_x.nprocs ) * sizeof(int));
    // disps_gatherv_x = malloc((comm_1d_x.nprocs ) * sizeof(int));
    // cnt_gatherv_y = malloc((comm_1d_y.nprocs ) * sizeof(int));
    // disps_gatherv_y = malloc((comm_1d_y.nprocs ) * sizeof(int));
    // cnt_gatherv_z = malloc((comm_1d_z.nprocs ) * sizeof(int));
    // disps_gatherv_z = malloc((comm_1d_z.nprocs ) * sizeof(int));
    cnt_gatherv_x = (int*)malloc(comm_1d_x.nprocs*sizeof(int));
    disps_gatherv_x = (int*)malloc(comm_1d_x.nprocs*sizeof(int));
    cnt_gatherv_y = (int*)malloc(comm_1d_y.nprocs*sizeof(int));
    disps_gatherv_y = (int*)malloc(comm_1d_y.nprocs*sizeof(int));
    cnt_gatherv_z = (int*)malloc(comm_1d_z.nprocs*sizeof(int));
    disps_gatherv_z = (int*)malloc(comm_1d_z.nprocs*sizeof(int));


    // timer_stamp0(STAMP_AGG);
    if (lv_aggregation_x != 0)
    {
        nx_aggr = mg_sdm[lv_aggregation_x].nx;
        if (nx_aggr % comm_1d_x.nprocs != 0) {
            if (myrank == 0) fprintf(stderr, "[MG] invalid x aggregation division\n");
            MPI_Abort(comm_1d_x.mpi_comm, EXIT_FAILURE);
        }
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
        MPI_Type_free(&ddtype_temp1);
    
        starts[0] = 1;
        starts[1] = 1;
        starts[2] = 1;
        MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &ddtype_temp1);
        MPI_Type_size(MPI_DOUBLE, &r8size);
        lb = 0;
        extent = (MPI_Aint) r8size;
        MPI_Type_create_resized(ddtype_temp1, lb, extent, &ddtype_gatherv_x);
        MPI_Type_commit(&ddtype_gatherv_x);
        MPI_Type_free(&ddtype_temp1);

        for (i = 0; i < comm_1d_x.nprocs; i++) {
            cnt_gatherv_x[i] = 1;
            disps_gatherv_x[i] = nx * i * (ny+2) * (nz+2);    
        }
    }

    if (lv_aggregation_y != 0)
    {
        nx_aggr = mg_sdm[lv_aggregation_y].nx;
        ny_aggr = mg_sdm[lv_aggregation_y].ny;
        if (ny_aggr % comm_1d_y.nprocs != 0) {
            if (myrank == 0) fprintf(stderr, "[MG] invalid y aggregation division\n");
            MPI_Abort(comm_1d_y.mpi_comm, EXIT_FAILURE);
        }
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
        MPI_Type_free(&ddtype_temp1);
    
        starts[0] = 1;
        starts[1] = 1;
        starts[2] = 1;
        MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &ddtype_temp1);
        MPI_Type_size(MPI_DOUBLE, &r8size);
        lb = 0;

        extent = (MPI_Aint) r8size;
        MPI_Type_create_resized(ddtype_temp1, lb, extent, &ddtype_gatherv_y);
        MPI_Type_commit(&ddtype_gatherv_y);
        MPI_Type_free(&ddtype_temp1);

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
        if (nz_aggr % comm_1d_z.nprocs != 0) {
            if (myrank == 0) fprintf(stderr, "[MG] invalid z aggregation division\n");
            MPI_Abort(comm_1d_z.mpi_comm, EXIT_FAILURE);
        }
    
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
        MPI_Type_free(&ddtype_temp1);
    
        starts[0] = 1;
        starts[1] = 1;
        starts[2] = 1;
        MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &ddtype_temp1);
        MPI_Type_size(MPI_DOUBLE, &r8size);
        lb = 0;
        extent = (MPI_Aint) r8size;
        MPI_Type_create_resized(ddtype_temp1, lb, extent, &ddtype_gatherv_z);
        MPI_Type_commit(&ddtype_gatherv_z);
        MPI_Type_free(&ddtype_temp1);

        for (i = 0; i < comm_1d_z.nprocs; i++) {
            cnt_gatherv_z[i] = 1;
            disps_gatherv_z[i] = nz * i;    
        }
    }
    // timer_stamp(8,STAMP_AGG);

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
        // if(myrank==0) printf("myrank=%d, lv_aggr_max_ptr= %d, lv_aggr_med_ptr= %d, lv_aggr_min_ptr= %d, max=%c, med=%c, min=%c\n",myrank, *lv_aggr_max_ptr, *lv_aggr_med_ptr, *lv_aggr_min_ptr, max, med, min );
    }
    // if(myrank==0) printf("level_case=%d\n", level_case);

    nx = sdm->nx;
    ny = sdm->ny;
    nz = sdm->nz;
    size_t size = (size_t)(nx + 2) * (size_t)(ny + 2) * (size_t)(nz + 2);
    double *rsd = NULL;
    CUDA_CHECK(cudaMalloc((void **)&rsd, size * sizeof(double)));
    CUDA_CHECK(cudaMemset(rsd, 0, size * sizeof(double)));

    switch (level_case) {
        case 1:
            // 1 aggregation level

            // timer_stamp0(STAMP_residual);
            multigrid_residual_gpu(rsd, a_poisson->coeff, sol, rhs, sdm, sdm->is_aggregated);
            // timer_stamp(18,STAMP_residual);
            vv_dot_3d_matrix_gpu(&res0tol, rsd, rsd, nx, ny, nz, sdm->is_aggregated);

            for(cyc = 1; cyc <= n_vcycles; cyc++)
            {
                // timer_stamp0(STAMP_LEVEL);
                // timer_stamp0(STAMP_smooth);
                rbgs_iterator_poisson_matrix_gpu(sol, a_poisson->coeff, rhs, sdm, maxiteration, omega_sor, sdm->is_aggregated);
                // timer_stamp(15,STAMP_smooth);
                // timer_stamp0(STAMP_residual);
                multigrid_residual_gpu(rsd, a_poisson->coeff, sol, rhs, sdm, sdm->is_aggregated);
                // timer_stamp(18,STAMP_residual);
                // timer_stamp(5,STAMP_LEVEL);
                // timer_stamp0(STAMP_restriction);
                multigrid_restriction_gpu(mg_sdm[1].b, rsd, &mg_sdm[1], sdm, 0);
                // timer_stamp(16,STAMP_restriction);
                // if(myrank==0) printf("[MG] Restriction from level 0 to level 1\n");

                for(l=1; l<=*lv_aggr_max_ptr-1; l++)
                {
                    size_t level_size = (size_t)(mg_sdm[l].nx + 2) *
                                        (size_t)(mg_sdm[l].ny + 2) *
                                        (size_t)(mg_sdm[l].nz + 2);
                    CUDA_CHECK(cudaMemset(mg_sdm[l].x, 0, level_size * sizeof(double)));
                    // timer_stamp0(STAMP_LEVEL);
                    // timer_stamp0(STAMP_smooth);
                    rbgs_iterator_poisson_matrix_gpu(mg_sdm[l].x, mg_a_poisson[l].coeff, mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                    // timer_stamp(15,STAMP_smooth);
                    // timer_stamp0(STAMP_residual);
                    multigrid_residual_gpu(mg_sdm[l].r, mg_a_poisson[l].coeff, mg_sdm[l].x, mg_sdm[l].b, &mg_sdm[l], mg_sdm[l].is_aggregated);
                    // timer_stamp(18,STAMP_residual);
                    // timer_stamp(6,STAMP_LEVEL);
                    // timer_stamp0(STAMP_restriction);
                    multigrid_restriction_gpu(mg_sdm[l+1].b, mg_sdm[l].r, &mg_sdm[l+1], &mg_sdm[l], l);
                    // timer_stamp(16,STAMP_restriction);
                    // if(myrank==0) printf("[MG] Restriction from level %d to level %d\n", l, l+1);
                }
    
l = *lv_aggr_max_ptr;
// timer_stamp0(STAMP_AGG);
if (comm_max_ptr->myrank == 0) {
    // timer_stamp0(STAMP_COMM_GATHER);
    CUDA_CHECK(cudaDeviceSynchronize());

    MPI_Gatherv(MPI_IN_PLACE, 0, *ddtype_send_max_ptr, mg_sdm[l].b, cnt_gatherv_max_ptr, disps_gatherv_max_ptr, *ddtype_gatherv_max_ptr, 0, comm_max_ptr->mpi_comm);
    // timer_stamp(13,STAMP_COMM_GATHER);
} else {
    // timer_stamp0(STAMP_COMM_GATHER);
    CUDA_CHECK(cudaDeviceSynchronize());

    MPI_Gatherv(mg_sdm[l].b, 1, *ddtype_send_max_ptr, mg_sdm[l].b, cnt_gatherv_max_ptr, disps_gatherv_max_ptr, *ddtype_gatherv_max_ptr, 0, comm_max_ptr->mpi_comm);
    // timer_stamp(13,STAMP_COMM_GATHER);
}
// timer_stamp(8,STAMP_AGG);

                if ((*comm_max_ptr).myrank == 0) {
                    // Restriction phase
                    for (int l = *lv_aggr_max_ptr; l < n_levels; l++) {
                    size_t level_size = (size_t)(mg_sdm[l].nx + 2) *
                                        (size_t)(mg_sdm[l].ny + 2) *
                                        (size_t)(mg_sdm[l].nz + 2);
                        CUDA_CHECK(cudaMemset(mg_sdm[l].x, 0, level_size * sizeof(double)));
                        // timer_stamp0(STAMP_LEVEL);
                        // timer_stamp0(STAMP_smooth);
                        rbgs_iterator_poisson_matrix_gpu(mg_sdm[l].x, mg_a_poisson[l].coeff, mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                        // timer_stamp(15,STAMP_smooth);
                        // timer_stamp0(STAMP_residual);
                        multigrid_residual_gpu(mg_sdm[l].r, mg_a_poisson[l].coeff, mg_sdm[l].x, mg_sdm[l].b, &mg_sdm[l], mg_sdm[l].is_aggregated);
                        // timer_stamp(18,STAMP_residual);
                        // timer_stamp(6,STAMP_LEVEL);
                        // timer_stamp0(STAMP_restriction);
                        multigrid_restriction_gpu(mg_sdm[l+1].b, mg_sdm[l].r, &mg_sdm[l+1], &mg_sdm[l], l);
                        // timer_stamp(16,STAMP_restriction);
                        // printf("[MG] Restriction from level %d to level %d\n", l, l+1);
                    }

                    // Solve at coarsest level
                    // timer_stamp0(STAMP_LEVEL);
                    // timer_stamp0(STAMP_smooth);
                    multigrid_solve_coarset_level_gpu(mg_sdm[n_levels].x, &mg_a_poisson[n_levels], mg_sdm[n_levels].b, &mg_sdm[n_levels], 1000, tolerance, omega_sor, mg_sdm[n_levels].is_aggregated);
                    // timer_stamp(15,STAMP_smooth);
                    // timer_stamp0(STAMP_residual);
                    multigrid_residual_gpu(mg_sdm[n_levels].r, mg_a_poisson[n_levels].coeff, mg_sdm[n_levels].x, mg_sdm[n_levels].b, &mg_sdm[n_levels], mg_sdm[n_levels].is_aggregated);
                    // timer_stamp(18,STAMP_residual);
                    // timer_stamp(7,STAMP_LEVEL);
                    ni = (mg_sdm[n_levels].ny+2) * (mg_sdm[n_levels].nz+2);
                    nj = (mg_sdm[n_levels].nz+2);
                    // adaptive_print_coarsest_gpu(cyc);

/* DEBUG_COARSEST host printer omitted: arrays reside on device. */

                    // Prolongation phase
                    for (l = n_levels-1; l >= *lv_aggr_max_ptr; l--) {
                        // timer_stamp0(STAMP_COMP);
                        // timer_stamp0(STAMP_prolongation);
                        multigrid_prolongation_linear_on_nonuniform_grid_gpu(mg_sdm[l].r, mg_sdm[l+1].x, &mg_sdm[l], &mg_sdm[l+1], l);
                        // timer_stamp(17,STAMP_prolongation);
                        // timer_stamp0(STAMP_LEVEL);
                        // x = x + r
                        nx = mg_sdm[l].nx;
                        ny = mg_sdm[l].ny;
                        nz = mg_sdm[l].nz;
                        adaptive_vector_add_gpu(mg_sdm[l].x, mg_sdm[l].r, (size_t)(nx+2)*(ny+2)*(nz+2));
                        // timer_stamp(10,STAMP_COMP);

                        // timer_stamp0(STAMP_smooth);
                        rbgs_iterator_poisson_matrix_gpu(mg_sdm[l].x, mg_a_poisson[l].coeff, mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                        // timer_stamp(15,STAMP_smooth);
                        // timer_stamp(6,STAMP_LEVEL);
                    }
                }

l = *lv_aggr_max_ptr;

// timer_stamp0(STAMP_AGG);
// #ifdef MPI_IN_PLACE
if (comm_max_ptr->myrank == 0) {
    // timer_stamp0(STAMP_COMM_SCATTER);
    CUDA_CHECK(cudaDeviceSynchronize());

    MPI_Scatterv(mg_sdm[l].x, cnt_gatherv_max_ptr, disps_gatherv_max_ptr, *ddtype_gatherv_max_ptr, MPI_IN_PLACE, 0, *ddtype_send_max_ptr, 0, comm_max_ptr->mpi_comm);
    // timer_stamp(14,STAMP_COMM_SCATTER);
} else {
    // timer_stamp0(STAMP_COMM_SCATTER);
    CUDA_CHECK(cudaDeviceSynchronize());

    MPI_Scatterv(mg_sdm[l].x, cnt_gatherv_max_ptr, disps_gatherv_max_ptr, *ddtype_gatherv_max_ptr, mg_sdm[l].x, 1, *ddtype_send_max_ptr, 0, comm_max_ptr->mpi_comm);
    // timer_stamp(14,STAMP_COMM_SCATTER);
}
// timer_stamp(8,STAMP_AGG);


                // Loop for prolongation from aggregation level down to level 1
                for (int l = *lv_aggr_max_ptr-1; l >= 1; l--) {
                    // timer_stamp0(STAMP_COMM_NEIGHBOR);
                    geometry_halocell_update_selectively_gpu(mg_sdm[l+1].x, &mg_sdm[l+1], mg_sdm[l+1].is_aggregated);
                    // timer_stamp(11,STAMP_COMM_NEIGHBOR);

                    // timer_stamp0(STAMP_COMP);
                    // timer_stamp0(STAMP_prolongation);
                    multigrid_prolongation_linear_on_nonuniform_grid_gpu(mg_sdm[l].r, mg_sdm[l+1].x,
                                                     &mg_sdm[l], &mg_sdm[l+1], l);
                    // timer_stamp(17,STAMP_prolongation);
                    // timer_stamp0(STAMP_LEVEL);
                    nx = mg_sdm[l].nx;
                    ny = mg_sdm[l].ny;
                    nz = mg_sdm[l].nz;
                    size_t size = (size_t)(nx + 2) * (size_t)(ny + 2) * (size_t)(nz + 2);
                    adaptive_vector_add_gpu(mg_sdm[l].x, mg_sdm[l].r, (size_t)size);
                    // timer_stamp(10,STAMP_COMP);
                    // timer_stamp0(STAMP_smooth);
                    rbgs_iterator_poisson_matrix_gpu(mg_sdm[l].x, mg_a_poisson[l].coeff, mg_sdm[l].b,
                                 &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                    // timer_stamp(15,STAMP_smooth);
                    // timer_stamp(6,STAMP_LEVEL);

                }

                // Prolongation to actual solution
                // timer_stamp0(STAMP_COMM_NEIGHBOR);
                geometry_halocell_update_selectively_gpu(mg_sdm[1].x, &mg_sdm[1], mg_sdm[1].is_aggregated);
                // timer_stamp(11,STAMP_COMM_NEIGHBOR);

                // timer_stamp0(STAMP_COMP);
                // timer_stamp0(STAMP_prolongation);
                multigrid_prolongation_linear_on_nonuniform_grid_gpu(rsd, mg_sdm[1].x, sdm, &mg_sdm[1], 0);
                // timer_stamp(17,STAMP_prolongation);
                // timer_stamp0(STAMP_LEVEL);
                nx = sdm->nx;
                ny = sdm->ny;
                nz = sdm->nz;
                    size_t size = (size_t)(nx + 2) * (size_t)(ny + 2) * (size_t)(nz + 2);
                adaptive_vector_add_gpu(sol, rsd, (size_t)size);
                // timer_stamp(10,STAMP_COMP);
                // timer_stamp0(STAMP_smooth);


                rbgs_iterator_poisson_matrix_gpu(sol, a_poisson->coeff, rhs, sdm, maxiteration, omega_sor, sdm->is_aggregated);
                // timer_stamp(15,STAMP_smooth);
                // timer_stamp(5,STAMP_LEVEL);

                // timer_stamp0(STAMP_residual);
                multigrid_residual_gpu(rsd, a_poisson->coeff, sol, rhs, sdm, sdm->is_aggregated);
                // timer_stamp(18,STAMP_residual);

                vv_dot_3d_matrix_gpu(&rsd_val, rsd, rsd, nx, ny, nz, sdm->is_aggregated);

                // if (myrank == 0) {
                //     printf("[MG] cycle = %d, Error: %e %e %e\n", cyc, sqrt(rsd_val), sqrt(res0tol), sqrt(rsd_val/res0tol));
                // }

                if (sqrt(rsd_val/res0tol) < tolerance) break;


            }

            // if(myrank==0) printf("[MG] Total %d V-cycles end\n", cyc);

            CUDA_CHECK(cudaFree(rsd));


            break;

        case 2:
            // 2 aggregation level
            // timer_stamp0(STAMP_residual);
            multigrid_residual_gpu(rsd, a_poisson->coeff, sol, rhs, sdm, sdm->is_aggregated);
            // timer_stamp(18,STAMP_residual);
            vv_dot_3d_matrix_gpu(&res0tol, rsd, rsd, nx, ny, nz, sdm->is_aggregated);

            for(cyc = 1; cyc <= n_vcycles; cyc++)
            {
                // timer_stamp0(STAMP_LEVEL);
                // timer_stamp0(STAMP_smooth);
                rbgs_iterator_poisson_matrix_gpu(sol, a_poisson->coeff, rhs, sdm, maxiteration, omega_sor, sdm->is_aggregated);
                // timer_stamp(15,STAMP_smooth);
                // timer_stamp0(STAMP_residual);
                multigrid_residual_gpu(rsd, a_poisson->coeff, sol, rhs, sdm, sdm->is_aggregated);
                // timer_stamp(18,STAMP_residual);
                // timer_stamp(5,STAMP_LEVEL);

                // timer_stamp0(STAMP_restriction);
                multigrid_restriction_gpu(mg_sdm[1].b, rsd, &mg_sdm[1], sdm, 0);
                // timer_stamp(16,STAMP_restriction);
                // if(myrank==0) printf("[MG] Restriction from level 0 to level 1\n");

                for(l=1; l<=*lv_aggr_med_ptr-1; l++)
                {
                    size_t level_size = (size_t)(mg_sdm[l].nx + 2) *
                                        (size_t)(mg_sdm[l].ny + 2) *
                                        (size_t)(mg_sdm[l].nz + 2);
                    CUDA_CHECK(cudaMemset(mg_sdm[l].x, 0, level_size * sizeof(double)));
                    // timer_stamp0(STAMP_LEVEL);
                    // timer_stamp0(STAMP_smooth);
                    rbgs_iterator_poisson_matrix_gpu(mg_sdm[l].x, mg_a_poisson[l].coeff, mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                    // timer_stamp(15,STAMP_smooth);
                    // timer_stamp0(STAMP_residual);
                    multigrid_residual_gpu(mg_sdm[l].r, mg_a_poisson[l].coeff, mg_sdm[l].x, mg_sdm[l].b, &mg_sdm[l], mg_sdm[l].is_aggregated);
                    // timer_stamp(18,STAMP_residual);
                    // timer_stamp(6,STAMP_LEVEL);
                    
                    // timer_stamp0(STAMP_restriction);
                    multigrid_restriction_gpu(mg_sdm[l+1].b, mg_sdm[l].r, &mg_sdm[l+1], &mg_sdm[l], l);
                    // timer_stamp(16,STAMP_restriction);
                    // if(myrank==0) printf("[MG] Restriction from level %d to level %d\n", l, l+1);
                }
    
l = *lv_aggr_med_ptr;
// timer_stamp0(STAMP_AGG);
if (comm_med_ptr->myrank == 0) {
    // timer_stamp0(STAMP_COMM_GATHER);
    CUDA_CHECK(cudaDeviceSynchronize());

    MPI_Gatherv(MPI_IN_PLACE, 0, *ddtype_send_med_ptr, mg_sdm[l].b, cnt_gatherv_med_ptr, disps_gatherv_med_ptr, *ddtype_gatherv_med_ptr, 0, comm_med_ptr->mpi_comm);
    // timer_stamp(13,STAMP_COMM_GATHER);
} else {
    // timer_stamp0(STAMP_COMM_GATHER);
    CUDA_CHECK(cudaDeviceSynchronize());

    MPI_Gatherv(mg_sdm[l].b, 1, *ddtype_send_med_ptr, mg_sdm[l].b, cnt_gatherv_med_ptr, disps_gatherv_med_ptr, *ddtype_gatherv_med_ptr, 0, comm_med_ptr->mpi_comm);
    // timer_stamp(13,STAMP_COMM_GATHER);
}
// timer_stamp(8,STAMP_AGG);

                if ((*comm_med_ptr).myrank == 0) {
                    // Restriction phase
                    for (int l = *lv_aggr_med_ptr; l < *lv_aggr_max_ptr; l++) {
                    size_t level_size = (size_t)(mg_sdm[l].nx + 2) *
                                        (size_t)(mg_sdm[l].ny + 2) *
                                        (size_t)(mg_sdm[l].nz + 2);
                        CUDA_CHECK(cudaMemset(mg_sdm[l].x, 0, level_size * sizeof(double)));
                        // timer_stamp0(STAMP_LEVEL);
                        // timer_stamp0(STAMP_smooth);
                        rbgs_iterator_poisson_matrix_gpu(mg_sdm[l].x, mg_a_poisson[l].coeff, mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                        // timer_stamp(15,STAMP_smooth);
                        // timer_stamp0(STAMP_residual);
                        multigrid_residual_gpu(mg_sdm[l].r, mg_a_poisson[l].coeff, mg_sdm[l].x, mg_sdm[l].b, &mg_sdm[l], mg_sdm[l].is_aggregated);
                        // timer_stamp(18,STAMP_residual);
                        // timer_stamp(6,STAMP_LEVEL);

                        // timer_stamp0(STAMP_restriction);
                        multigrid_restriction_gpu(mg_sdm[l+1].b, mg_sdm[l].r, &mg_sdm[l+1], &mg_sdm[l], l);
                        // timer_stamp(16,STAMP_restriction);
                        // printf("[MG] Restriction from level %d to level %d\n", l, l+1);
                    }

l = *lv_aggr_max_ptr;
// timer_stamp0(STAMP_AGG);
if (comm_max_ptr->myrank == 0) {
    // timer_stamp0(STAMP_COMM_GATHER);
    CUDA_CHECK(cudaDeviceSynchronize());

    MPI_Gatherv(MPI_IN_PLACE, 0, *ddtype_send_max_ptr, mg_sdm[l].b, cnt_gatherv_max_ptr, disps_gatherv_max_ptr, *ddtype_gatherv_max_ptr, 0, comm_max_ptr->mpi_comm);
    // timer_stamp(13,STAMP_COMM_GATHER);
} else {
    // timer_stamp0(STAMP_COMM_GATHER);
    CUDA_CHECK(cudaDeviceSynchronize());

    MPI_Gatherv(mg_sdm[l].b, 1, *ddtype_send_max_ptr, mg_sdm[l].b, cnt_gatherv_max_ptr, disps_gatherv_max_ptr, *ddtype_gatherv_max_ptr, 0, comm_max_ptr->mpi_comm);
    // timer_stamp(13,STAMP_COMM_GATHER);
}
// timer_stamp(8,STAMP_AGG);

                if ((*comm_max_ptr).myrank == 0) {
                    // Restriction phase
                    for (int l = *lv_aggr_max_ptr; l < n_levels; l++) {
                    size_t level_size = (size_t)(mg_sdm[l].nx + 2) *
                                        (size_t)(mg_sdm[l].ny + 2) *
                                        (size_t)(mg_sdm[l].nz + 2);
                        CUDA_CHECK(cudaMemset(mg_sdm[l].x, 0, level_size * sizeof(double)));
                        // timer_stamp0(STAMP_LEVEL);
                        // timer_stamp0(STAMP_smooth);
                        rbgs_iterator_poisson_matrix_gpu(mg_sdm[l].x, mg_a_poisson[l].coeff, mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                        // timer_stamp(15,STAMP_smooth);
                        // timer_stamp0(STAMP_residual);
                        multigrid_residual_gpu(mg_sdm[l].r, mg_a_poisson[l].coeff, mg_sdm[l].x, mg_sdm[l].b, &mg_sdm[l], mg_sdm[l].is_aggregated);
                        // timer_stamp(18,STAMP_residual);
                        // timer_stamp(6,STAMP_LEVEL);
                        
                        // timer_stamp0(STAMP_restriction);
                        multigrid_restriction_gpu(mg_sdm[l+1].b, mg_sdm[l].r, &mg_sdm[l+1], &mg_sdm[l], l);
                        // timer_stamp(16,STAMP_restriction);
                        // printf("[MG] Restriction from level %d to level %d\n", l, l+1);
                    }

                    // Solve at coarsest level
                    // timer_stamp0(STAMP_LEVEL);
                    // timer_stamp0(STAMP_smooth);
                    multigrid_solve_coarset_level_gpu(mg_sdm[n_levels].x, &mg_a_poisson[n_levels], mg_sdm[n_levels].b, &mg_sdm[n_levels], 1000, tolerance, omega_sor, mg_sdm[n_levels].is_aggregated);
                    // timer_stamp(15,STAMP_smooth);
                    // timer_stamp0(STAMP_residual);
                    multigrid_residual_gpu(mg_sdm[n_levels].r, mg_a_poisson[n_levels].coeff, mg_sdm[n_levels].x, mg_sdm[n_levels].b, &mg_sdm[n_levels], mg_sdm[n_levels].is_aggregated);
                    // timer_stamp(18,STAMP_residual);
                    // timer_stamp(7,STAMP_LEVEL);

                    ni = (mg_sdm[n_levels].ny+2) * (mg_sdm[n_levels].nz+2);
                    nj = (mg_sdm[n_levels].nz+2);
                    // adaptive_print_coarsest_gpu(cyc);

/* DEBUG_COARSEST host printer omitted: arrays reside on device. */

                    // Prolongation phase
                    for (l = n_levels-1; l >= *lv_aggr_max_ptr; l--) {
                        // timer_stamp0(STAMP_COMP);
                        // timer_stamp0(STAMP_prolongation);
                        multigrid_prolongation_linear_on_nonuniform_grid_gpu(mg_sdm[l].r, mg_sdm[l+1].x, &mg_sdm[l], &mg_sdm[l+1], l);
                        // timer_stamp(17,STAMP_prolongation);
                        // timer_stamp0(STAMP_LEVEL);
                        // x = x + r
                        nx = mg_sdm[l].nx;
                        ny = mg_sdm[l].ny;
                        nz = mg_sdm[l].nz;
                        adaptive_vector_add_gpu(mg_sdm[l].x, mg_sdm[l].r, (size_t)(nx+2)*(ny+2)*(nz+2));
                        // timer_stamp(10,STAMP_COMP);

                    // timer_stamp0(STAMP_smooth);
                    rbgs_iterator_poisson_matrix_gpu(mg_sdm[l].x, mg_a_poisson[l].coeff, mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                    // timer_stamp(15,STAMP_smooth);
                    // timer_stamp(6,STAMP_LEVEL);
                    }
                }

l = *lv_aggr_max_ptr;
// timer_stamp0(STAMP_AGG);
// #ifdef MPI_IN_PLACE
if (comm_max_ptr->myrank == 0) {
    // timer_stamp0(STAMP_COMM_SCATTER);
    CUDA_CHECK(cudaDeviceSynchronize());

    MPI_Scatterv(mg_sdm[l].x, cnt_gatherv_max_ptr, disps_gatherv_max_ptr, *ddtype_gatherv_max_ptr, MPI_IN_PLACE, 0, *ddtype_send_max_ptr, 0, comm_max_ptr->mpi_comm);

    // timer_stamp(14,STAMP_COMM_SCATTER);
} else {
    // timer_stamp0(STAMP_COMM_SCATTER);
    CUDA_CHECK(cudaDeviceSynchronize());

    MPI_Scatterv(mg_sdm[l].x, cnt_gatherv_max_ptr, disps_gatherv_max_ptr, *ddtype_gatherv_max_ptr, mg_sdm[l].x, 1, *ddtype_send_max_ptr, 0, comm_max_ptr->mpi_comm);
    // timer_stamp(14,STAMP_COMM_SCATTER);
}
// timer_stamp(8,STAMP_AGG);


                // Loop for prolongation from aggregation level down to level 1
                for (int l = *lv_aggr_max_ptr-1; l >= *lv_aggr_med_ptr; l--) {
                    // timer_stamp0(STAMP_COMM_NEIGHBOR);
                    geometry_halocell_update_selectively_gpu(mg_sdm[l+1].x, &mg_sdm[l+1], mg_sdm[l+1].is_aggregated);
                    // timer_stamp(11,STAMP_COMM_NEIGHBOR);

                    // timer_stamp0(STAMP_COMP);
                    // timer_stamp0(STAMP_prolongation);
                    multigrid_prolongation_linear_on_nonuniform_grid_gpu(mg_sdm[l].r, mg_sdm[l+1].x,
                                                     &mg_sdm[l], &mg_sdm[l+1], l);
                    // timer_stamp(17,STAMP_prolongation);
                    // timer_stamp0(STAMP_LEVEL);
                    nx = mg_sdm[l].nx;
                    ny = mg_sdm[l].ny;
                    nz = mg_sdm[l].nz;
                    size_t size = (size_t)(nx + 2) * (size_t)(ny + 2) * (size_t)(nz + 2);
                    adaptive_vector_add_gpu(mg_sdm[l].x, mg_sdm[l].r, (size_t)size);
                    // timer_stamp(10,STAMP_COMP);
                    // timer_stamp0(STAMP_smooth);
                    rbgs_iterator_poisson_matrix_gpu(mg_sdm[l].x, mg_a_poisson[l].coeff, mg_sdm[l].b,
                                 &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                    // timer_stamp(15,STAMP_smooth);
                    // timer_stamp(6,STAMP_LEVEL);
                }
            }


l = *lv_aggr_med_ptr;
// timer_stamp0(STAMP_AGG);
// #ifdef MPI_IN_PLACE
if (comm_med_ptr->myrank == 0) {
    // timer_stamp0(STAMP_COMM_SCATTER);
    CUDA_CHECK(cudaDeviceSynchronize());

    MPI_Scatterv(mg_sdm[l].x, cnt_gatherv_med_ptr, disps_gatherv_med_ptr, *ddtype_gatherv_med_ptr, MPI_IN_PLACE, 0, *ddtype_send_med_ptr, 0, comm_med_ptr->mpi_comm);
    // timer_stamp(14,STAMP_COMM_SCATTER);
} else {
    // timer_stamp0(STAMP_COMM_SCATTER);
    CUDA_CHECK(cudaDeviceSynchronize());

    MPI_Scatterv(mg_sdm[l].x, cnt_gatherv_med_ptr, disps_gatherv_med_ptr, *ddtype_gatherv_med_ptr, mg_sdm[l].x, 1, *ddtype_send_med_ptr, 0, comm_med_ptr->mpi_comm);
    // timer_stamp(14,STAMP_COMM_SCATTER);
}
// timer_stamp(8,STAMP_AGG);


                // Loop for prolongation from aggregation level down to level 1
                for (int l = *lv_aggr_med_ptr-1; l >= 1; l--) {
                    // timer_stamp0(STAMP_COMM_NEIGHBOR);
                    geometry_halocell_update_selectively_gpu(mg_sdm[l+1].x, &mg_sdm[l+1], mg_sdm[l+1].is_aggregated);
                    // timer_stamp(11,STAMP_COMM_NEIGHBOR);

                    // timer_stamp0(STAMP_COMP);
                    // timer_stamp0(STAMP_prolongation);
                    multigrid_prolongation_linear_on_nonuniform_grid_gpu(mg_sdm[l].r, mg_sdm[l+1].x,
                                                     &mg_sdm[l], &mg_sdm[l+1], l);

                    // timer_stamp(17,STAMP_prolongation);
                    
                    // timer_stamp0(STAMP_LEVEL);
                    nx = mg_sdm[l].nx;
                    ny = mg_sdm[l].ny;
                    nz = mg_sdm[l].nz;
                    size_t size = (size_t)(nx + 2) * (size_t)(ny + 2) * (size_t)(nz + 2);
                    adaptive_vector_add_gpu(mg_sdm[l].x, mg_sdm[l].r, (size_t)size);
                    // timer_stamp(10,STAMP_COMP);
                    // timer_stamp0(STAMP_smooth);
                    rbgs_iterator_poisson_matrix_gpu(mg_sdm[l].x, mg_a_poisson[l].coeff, mg_sdm[l].b,
                                 &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                    // timer_stamp(15,STAMP_smooth);
                    // timer_stamp(6,STAMP_LEVEL);
                }

                // Prolongation to actual solution
                // timer_stamp0(STAMP_COMM_NEIGHBOR);
                geometry_halocell_update_selectively_gpu(mg_sdm[1].x, &mg_sdm[1], mg_sdm[1].is_aggregated);
                // timer_stamp(11,STAMP_COMM_NEIGHBOR);

                // timer_stamp0(STAMP_COMP);
                // timer_stamp0(STAMP_prolongation);
                multigrid_prolongation_linear_on_nonuniform_grid_gpu(rsd, mg_sdm[1].x, sdm, &mg_sdm[1], 0);
                // timer_stamp(17,STAMP_prolongation);
                // timer_stamp0(STAMP_LEVEL);
                nx = sdm->nx;
                ny = sdm->ny;
                nz = sdm->nz;
                    size_t size = (size_t)(nx + 2) * (size_t)(ny + 2) * (size_t)(nz + 2);
                adaptive_vector_add_gpu(sol, rsd, (size_t)size);
                // timer_stamp(10,STAMP_COMP);
                // timer_stamp0(STAMP_smooth);
                rbgs_iterator_poisson_matrix_gpu(sol, a_poisson->coeff, rhs, sdm, maxiteration, omega_sor, sdm->is_aggregated);
                // timer_stamp(15,STAMP_smooth);
                // timer_stamp(5,STAMP_LEVEL);

                // timer_stamp0(STAMP_residual);
                multigrid_residual_gpu(rsd, a_poisson->coeff, sol, rhs, sdm, sdm->is_aggregated);
                // timer_stamp(18,STAMP_residual);

                vv_dot_3d_matrix_gpu(&rsd_val, rsd, rsd, nx, ny, nz, sdm->is_aggregated);

                // if (myrank == 0) {
                //     printf("[MG] cycle = %d, Error: %e %e %e\n", cyc, sqrt(rsd_val), sqrt(res0tol), sqrt(rsd_val/res0tol));
                // }

                if (sqrt(rsd_val/res0tol) < tolerance) break;


            }

            // if(myrank==0) printf("[MG] Total %d V-cycles end\n", cyc);

            CUDA_CHECK(cudaFree(rsd));


            break;

        case 3:
            // 3 aggregation level

            // timer_stamp0(STAMP_residual);
            multigrid_residual_gpu(rsd, a_poisson->coeff, sol, rhs, sdm, sdm->is_aggregated);
            // timer_stamp(18,STAMP_residual);
            vv_dot_3d_matrix_gpu(&res0tol, rsd, rsd, nx, ny, nz, sdm->is_aggregated);

            for(cyc = 1; cyc <= n_vcycles; cyc++)
            {
                // timer_stamp0(STAMP_LEVEL);
                // timer_stamp0(STAMP_smooth);
                rbgs_iterator_poisson_matrix_gpu(sol, a_poisson->coeff, rhs, sdm, maxiteration, omega_sor, sdm->is_aggregated);
                // timer_stamp(15,STAMP_smooth);
                // timer_stamp0(STAMP_residual);
                multigrid_residual_gpu(rsd, a_poisson->coeff, sol, rhs, sdm, sdm->is_aggregated);
                // timer_stamp(18,STAMP_residual);
                // timer_stamp(5,STAMP_LEVEL);

                // timer_stamp0(STAMP_restriction);
                multigrid_restriction_gpu(mg_sdm[1].b, rsd, &mg_sdm[1], sdm, 0);
                // timer_stamp(16,STAMP_restriction);
                // if(myrank==0) printf("[MG] Restriction from level 0 to level 1\n");

                for(l=1; l<=*lv_aggr_min_ptr-1; l++)
                {
                    size_t level_size = (size_t)(mg_sdm[l].nx + 2) *
                                        (size_t)(mg_sdm[l].ny + 2) *
                                        (size_t)(mg_sdm[l].nz + 2);
                    CUDA_CHECK(cudaMemset(mg_sdm[l].x, 0, level_size * sizeof(double)));
                    // timer_stamp0(STAMP_LEVEL);
                    // timer_stamp0(STAMP_smooth);
                    rbgs_iterator_poisson_matrix_gpu(mg_sdm[l].x, mg_a_poisson[l].coeff, mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                    // timer_stamp(15,STAMP_smooth);
                    // timer_stamp0(STAMP_residual);
                    multigrid_residual_gpu(mg_sdm[l].r, mg_a_poisson[l].coeff, mg_sdm[l].x, mg_sdm[l].b, &mg_sdm[l], mg_sdm[l].is_aggregated);
                    // timer_stamp(18,STAMP_residual);
                    // timer_stamp(6,STAMP_LEVEL);
                    
                    // timer_stamp0(STAMP_restriction);
                    multigrid_restriction_gpu(mg_sdm[l+1].b, mg_sdm[l].r, &mg_sdm[l+1], &mg_sdm[l], l);
                    // timer_stamp(16,STAMP_restriction);


                    // if(myrank==0) printf("[MG] Restriction from level %d to level %d\n", l, l+1);
                }

l = *lv_aggr_min_ptr;
// timer_stamp0(STAMP_AGG);
if (comm_min_ptr->myrank == 0) {
    // timer_stamp0(STAMP_COMM_GATHER);
    CUDA_CHECK(cudaDeviceSynchronize());

    MPI_Gatherv(MPI_IN_PLACE, 0, *ddtype_send_min_ptr, mg_sdm[l].b, cnt_gatherv_min_ptr, disps_gatherv_min_ptr, *ddtype_gatherv_min_ptr, 0, comm_min_ptr->mpi_comm);
    // timer_stamp(13,STAMP_COMM_GATHER);
} else {
    // timer_stamp0(STAMP_COMM_GATHER);
    CUDA_CHECK(cudaDeviceSynchronize());

    MPI_Gatherv(mg_sdm[l].b, 1, *ddtype_send_min_ptr, mg_sdm[l].b, cnt_gatherv_min_ptr, disps_gatherv_min_ptr, *ddtype_gatherv_min_ptr, 0, comm_min_ptr->mpi_comm);
    // timer_stamp(13,STAMP_COMM_GATHER);
}
// timer_stamp(8,STAMP_AGG);

                if ((*comm_min_ptr).myrank == 0) {

                    for (int l = *lv_aggr_min_ptr; l < *lv_aggr_med_ptr; l++) {
                    size_t level_size = (size_t)(mg_sdm[l].nx + 2) *
                                        (size_t)(mg_sdm[l].ny + 2) *
                                        (size_t)(mg_sdm[l].nz + 2);
                        CUDA_CHECK(cudaMemset(mg_sdm[l].x, 0, level_size * sizeof(double)));
                        // timer_stamp0(STAMP_LEVEL);
                        // timer_stamp0(STAMP_smooth);
                        rbgs_iterator_poisson_matrix_gpu(mg_sdm[l].x, mg_a_poisson[l].coeff, mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                        // timer_stamp(15,STAMP_smooth);
                        // timer_stamp0(STAMP_residual);
                        multigrid_residual_gpu(mg_sdm[l].r, mg_a_poisson[l].coeff, mg_sdm[l].x, mg_sdm[l].b, &mg_sdm[l], mg_sdm[l].is_aggregated);
                        // timer_stamp(18,STAMP_residual);
                        // timer_stamp(6,STAMP_LEVEL);
                        
                        // timer_stamp0(STAMP_restriction);
                        multigrid_restriction_gpu(mg_sdm[l+1].b, mg_sdm[l].r, &mg_sdm[l+1], &mg_sdm[l], l);
                        // timer_stamp(16,STAMP_restriction);
                        // printf("[MG] Restriction from level %d to level %d\n", l, l+1);
                    }
    
l = *lv_aggr_med_ptr;
// timer_stamp0(STAMP_AGG);
if (comm_med_ptr->myrank == 0) {
    // timer_stamp0(STAMP_COMM_GATHER);
    CUDA_CHECK(cudaDeviceSynchronize());

    MPI_Gatherv(MPI_IN_PLACE, 0, *ddtype_send_med_ptr, mg_sdm[l].b, cnt_gatherv_med_ptr, disps_gatherv_med_ptr, *ddtype_gatherv_med_ptr, 0, comm_med_ptr->mpi_comm);
    // timer_stamp(13,STAMP_COMM_GATHER);
} else {
    // timer_stamp0(STAMP_COMM_GATHER);
    CUDA_CHECK(cudaDeviceSynchronize());

    MPI_Gatherv(mg_sdm[l].b, 1, *ddtype_send_med_ptr, mg_sdm[l].b, cnt_gatherv_med_ptr, disps_gatherv_med_ptr, *ddtype_gatherv_med_ptr, 0, comm_med_ptr->mpi_comm);
    // timer_stamp(13,STAMP_COMM_GATHER);
}
// timer_stamp(8,STAMP_AGG);

                if ((*comm_med_ptr).myrank == 0) {
                    // Restriction phase
                    for (int l = *lv_aggr_med_ptr; l < *lv_aggr_max_ptr; l++) {
                    size_t level_size = (size_t)(mg_sdm[l].nx + 2) *
                                        (size_t)(mg_sdm[l].ny + 2) *
                                        (size_t)(mg_sdm[l].nz + 2);
                        CUDA_CHECK(cudaMemset(mg_sdm[l].x, 0, level_size * sizeof(double)));
                        // timer_stamp0(STAMP_LEVEL);
                        // timer_stamp0(STAMP_smooth);
                        rbgs_iterator_poisson_matrix_gpu(mg_sdm[l].x, mg_a_poisson[l].coeff, mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                        // timer_stamp(15,STAMP_smooth);
                        // timer_stamp0(STAMP_residual);
                        multigrid_residual_gpu(mg_sdm[l].r, mg_a_poisson[l].coeff, mg_sdm[l].x, mg_sdm[l].b, &mg_sdm[l], mg_sdm[l].is_aggregated);
                        // timer_stamp(18,STAMP_residual);
                        // timer_stamp(6,STAMP_LEVEL);

                        // timer_stamp0(STAMP_restriction);
                        multigrid_restriction_gpu(mg_sdm[l+1].b, mg_sdm[l].r, &mg_sdm[l+1], &mg_sdm[l], l);
                        // timer_stamp(16,STAMP_restriction);
                        // printf("[MG] Restriction from level %d to level %d\n", l, l+1);
                    }

l = *lv_aggr_max_ptr;
// timer_stamp0(STAMP_AGG);
if (comm_max_ptr->myrank == 0) {
    // timer_stamp0(STAMP_COMM_GATHER);

    CUDA_CHECK(cudaDeviceSynchronize());

    MPI_Gatherv(MPI_IN_PLACE, 0, *ddtype_send_max_ptr, mg_sdm[l].b, cnt_gatherv_max_ptr, disps_gatherv_max_ptr, *ddtype_gatherv_max_ptr, 0, comm_max_ptr->mpi_comm);
    // timer_stamp(13,STAMP_COMM_GATHER);
} else {
    // timer_stamp0(STAMP_COMM_GATHER);
    CUDA_CHECK(cudaDeviceSynchronize());

    MPI_Gatherv(mg_sdm[l].b, 1, *ddtype_send_max_ptr, mg_sdm[l].b, cnt_gatherv_max_ptr, disps_gatherv_max_ptr, *ddtype_gatherv_max_ptr, 0, comm_max_ptr->mpi_comm);
    // timer_stamp(13,STAMP_COMM_GATHER);
}
// timer_stamp(8,STAMP_AGG);

                if ((*comm_max_ptr).myrank == 0) {
                    // Restriction phase
                    for (int l = *lv_aggr_max_ptr; l < n_levels; l++) {
                    size_t level_size = (size_t)(mg_sdm[l].nx + 2) *
                                        (size_t)(mg_sdm[l].ny + 2) *
                                        (size_t)(mg_sdm[l].nz + 2);
                        CUDA_CHECK(cudaMemset(mg_sdm[l].x, 0, level_size * sizeof(double)));
                        // timer_stamp0(STAMP_LEVEL);
                        // timer_stamp0(STAMP_smooth);
                        rbgs_iterator_poisson_matrix_gpu(mg_sdm[l].x, mg_a_poisson[l].coeff, mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                        // timer_stamp(15,STAMP_smooth);
                        // timer_stamp0(STAMP_residual);
                        multigrid_residual_gpu(mg_sdm[l].r, mg_a_poisson[l].coeff, mg_sdm[l].x, mg_sdm[l].b, &mg_sdm[l], mg_sdm[l].is_aggregated);
                        // timer_stamp(18,STAMP_residual);
                        // timer_stamp(6,STAMP_LEVEL);
                        // timer_stamp0(STAMP_restriction);
                        multigrid_restriction_gpu(mg_sdm[l+1].b, mg_sdm[l].r, &mg_sdm[l+1], &mg_sdm[l], l);
                        // timer_stamp(16,STAMP_restriction);
                        // printf("[MG] Restriction from level %d to level %d\n", l, l+1);
                    }

                    // Solve at coarsest level
                    // timer_stamp0(STAMP_LEVEL);
                    // timer_stamp0(STAMP_smooth);
                    multigrid_solve_coarset_level_gpu(mg_sdm[n_levels].x, &mg_a_poisson[n_levels], mg_sdm[n_levels].b, &mg_sdm[n_levels], 1000, tolerance, omega_sor, mg_sdm[n_levels].is_aggregated);
                    // timer_stamp(15,STAMP_smooth);
                    // timer_stamp0(STAMP_residual);
                    multigrid_residual_gpu(mg_sdm[n_levels].r, mg_a_poisson[n_levels].coeff, mg_sdm[n_levels].x, mg_sdm[n_levels].b, &mg_sdm[n_levels], mg_sdm[n_levels].is_aggregated);
                    // timer_stamp(18,STAMP_residual);
                    // timer_stamp(7,STAMP_LEVEL);
                    ni = (mg_sdm[n_levels].ny+2) * (mg_sdm[n_levels].nz+2);
                    nj = (mg_sdm[n_levels].nz+2);
                    // adaptive_print_coarsest_gpu(cyc);

/* DEBUG_COARSEST host printer omitted: arrays reside on device. */

                    // Prolongation phase
                    for (l = n_levels-1; l >= *lv_aggr_max_ptr; l--) {
                        // timer_stamp0(STAMP_COMP);
                        // timer_stamp0(STAMP_prolongation);
                        multigrid_prolongation_linear_on_nonuniform_grid_gpu(mg_sdm[l].r, mg_sdm[l+1].x, &mg_sdm[l], &mg_sdm[l+1], l);
                        // timer_stamp(17,STAMP_prolongation);
                        // timer_stamp0(STAMP_LEVEL);
                        // x = x + r
                        nx = mg_sdm[l].nx;
                        ny = mg_sdm[l].ny;
                        nz = mg_sdm[l].nz;
                        adaptive_vector_add_gpu(mg_sdm[l].x, mg_sdm[l].r, (size_t)(nx+2)*(ny+2)*(nz+2));
                        // timer_stamp(10,STAMP_COMP);
                        // timer_stamp0(STAMP_smooth);
                        rbgs_iterator_poisson_matrix_gpu(mg_sdm[l].x, mg_a_poisson[l].coeff, mg_sdm[l].b, &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                        // timer_stamp(15,STAMP_smooth);
                        // timer_stamp(6,STAMP_LEVEL);
                    }
                }

l = *lv_aggr_max_ptr;
// timer_stamp0(STAMP_AGG);
// #ifdef MPI_IN_PLACE
if (comm_max_ptr->myrank == 0) {
    // timer_stamp0(STAMP_COMM_SCATTER);
    CUDA_CHECK(cudaDeviceSynchronize());

    MPI_Scatterv(mg_sdm[l].x, cnt_gatherv_max_ptr, disps_gatherv_max_ptr, *ddtype_gatherv_max_ptr, MPI_IN_PLACE, 0, *ddtype_send_max_ptr, 0, comm_max_ptr->mpi_comm);
    // timer_stamp(14,STAMP_COMM_SCATTER);
} else {
    // timer_stamp0(STAMP_COMM_SCATTER);
    CUDA_CHECK(cudaDeviceSynchronize());

    MPI_Scatterv(mg_sdm[l].x, cnt_gatherv_max_ptr, disps_gatherv_max_ptr, *ddtype_gatherv_max_ptr, mg_sdm[l].x, 1, *ddtype_send_max_ptr, 0, comm_max_ptr->mpi_comm);
    // timer_stamp(14,STAMP_COMM_SCATTER);
}
// timer_stamp(8,STAMP_AGG);

                // Loop for prolongation from aggregation level down to level 1
                for (int l = *lv_aggr_max_ptr-1; l >= *lv_aggr_med_ptr; l--) {
                    // timer_stamp0(STAMP_COMM_NEIGHBOR);
                    geometry_halocell_update_selectively_gpu(mg_sdm[l+1].x, &mg_sdm[l+1], mg_sdm[l+1].is_aggregated);
                    // timer_stamp(11,STAMP_COMM_NEIGHBOR);

                    // timer_stamp0(STAMP_COMP);
                    // timer_stamp0(STAMP_prolongation);
                    multigrid_prolongation_linear_on_nonuniform_grid_gpu(mg_sdm[l].r, mg_sdm[l+1].x,
                                                     &mg_sdm[l], &mg_sdm[l+1], l);
                    // timer_stamp(17,STAMP_prolongation);

                    // timer_stamp0(STAMP_LEVEL);
                    nx = mg_sdm[l].nx;
                    ny = mg_sdm[l].ny;
                    nz = mg_sdm[l].nz;
                    size_t size = (size_t)(nx + 2) * (size_t)(ny + 2) * (size_t)(nz + 2);
                    adaptive_vector_add_gpu(mg_sdm[l].x, mg_sdm[l].r, (size_t)size);
                    // timer_stamp(10,STAMP_COMP);
                    // timer_stamp0(STAMP_smooth);
                    rbgs_iterator_poisson_matrix_gpu(mg_sdm[l].x, mg_a_poisson[l].coeff, mg_sdm[l].b,
                                 &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                    // timer_stamp(15,STAMP_smooth);
                    // timer_stamp(6,STAMP_LEVEL);
                }
            }


l = *lv_aggr_med_ptr;
// timer_stamp0(STAMP_AGG);
// #ifdef MPI_IN_PLACE
if (comm_med_ptr->myrank == 0) {
    // timer_stamp0(STAMP_COMM_SCATTER);
    CUDA_CHECK(cudaDeviceSynchronize());

    MPI_Scatterv(mg_sdm[l].x, cnt_gatherv_med_ptr, disps_gatherv_med_ptr, *ddtype_gatherv_med_ptr, MPI_IN_PLACE, 0, *ddtype_send_med_ptr, 0, comm_med_ptr->mpi_comm);
    // timer_stamp(14,STAMP_COMM_SCATTER);
} else {
    // timer_stamp0(STAMP_COMM_SCATTER);
    CUDA_CHECK(cudaDeviceSynchronize());

    MPI_Scatterv(mg_sdm[l].x, cnt_gatherv_med_ptr, disps_gatherv_med_ptr, *ddtype_gatherv_med_ptr, mg_sdm[l].x, 1, *ddtype_send_med_ptr, 0, comm_med_ptr->mpi_comm);
    // timer_stamp(14,STAMP_COMM_SCATTER);
}
// timer_stamp(8,STAMP_AGG);

                for (int l = *lv_aggr_med_ptr-1; l >= *lv_aggr_min_ptr; l--) {
                    // timer_stamp0(STAMP_COMM_NEIGHBOR);
                    geometry_halocell_update_selectively_gpu(mg_sdm[l+1].x, &mg_sdm[l+1], mg_sdm[l+1].is_aggregated);
                    // timer_stamp(11,STAMP_COMM_NEIGHBOR);

                    // timer_stamp0(STAMP_COMP);
                    // timer_stamp0(STAMP_prolongation);
                    multigrid_prolongation_linear_on_nonuniform_grid_gpu(mg_sdm[l].r, mg_sdm[l+1].x,
                                                     &mg_sdm[l], &mg_sdm[l+1], l);
                    // timer_stamp(17,STAMP_prolongation);
                    // timer_stamp0(STAMP_LEVEL);
                    nx = mg_sdm[l].nx;
                    ny = mg_sdm[l].ny;
                    nz = mg_sdm[l].nz;
                    size_t size = (size_t)(nx + 2) * (size_t)(ny + 2) * (size_t)(nz + 2);
                    adaptive_vector_add_gpu(mg_sdm[l].x, mg_sdm[l].r, (size_t)size);
                    // timer_stamp(10,STAMP_COMP);
                    // timer_stamp0(STAMP_smooth);
                    rbgs_iterator_poisson_matrix_gpu(mg_sdm[l].x, mg_a_poisson[l].coeff, mg_sdm[l].b,
                                 &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                    // timer_stamp(15,STAMP_smooth);
                    // timer_stamp(6,STAMP_LEVEL);
                }
            }

l = *lv_aggr_min_ptr;
// timer_stamp0(STAMP_AGG);
// #ifdef MPI_IN_PLACE
if (comm_min_ptr->myrank == 0) {
    // timer_stamp0(STAMP_COMM_SCATTER);
    CUDA_CHECK(cudaDeviceSynchronize());

    MPI_Scatterv(mg_sdm[l].x, cnt_gatherv_min_ptr, disps_gatherv_min_ptr, *ddtype_gatherv_min_ptr, MPI_IN_PLACE, 0, *ddtype_send_min_ptr, 0, comm_min_ptr->mpi_comm);
    // timer_stamp(14,STAMP_COMM_SCATTER);
} else {
    // timer_stamp0(STAMP_COMM_SCATTER);
    CUDA_CHECK(cudaDeviceSynchronize());

    MPI_Scatterv(mg_sdm[l].x, cnt_gatherv_min_ptr, disps_gatherv_min_ptr, *ddtype_gatherv_min_ptr, mg_sdm[l].x, 1, *ddtype_send_min_ptr, 0, comm_min_ptr->mpi_comm);
    // timer_stamp(14,STAMP_COMM_SCATTER);
}
// timer_stamp(8,STAMP_AGG);

                // Loop for prolongation from aggregation level down to level 1
                for (int l = *lv_aggr_min_ptr-1; l >= 1; l--) {
                    // timer_stamp0(STAMP_COMM_NEIGHBOR);
                    geometry_halocell_update_selectively_gpu(mg_sdm[l+1].x, &mg_sdm[l+1], mg_sdm[l+1].is_aggregated);
                    // timer_stamp(11,STAMP_COMM_NEIGHBOR);

                    // timer_stamp0(STAMP_COMP);
                    // timer_stamp0(STAMP_prolongation);
                    multigrid_prolongation_linear_on_nonuniform_grid_gpu(mg_sdm[l].r, mg_sdm[l+1].x,
                                                     &mg_sdm[l], &mg_sdm[l+1], l);
                    // timer_stamp(17,STAMP_prolongation);

                    // timer_stamp0(STAMP_LEVEL);
                    nx = mg_sdm[l].nx;
                    ny = mg_sdm[l].ny;
                    nz = mg_sdm[l].nz;
                    size_t size = (size_t)(nx + 2) * (size_t)(ny + 2) * (size_t)(nz + 2);
                    adaptive_vector_add_gpu(mg_sdm[l].x, mg_sdm[l].r, (size_t)size);
                    // timer_stamp(10,STAMP_COMP);
                    // timer_stamp0(STAMP_smooth);
                    rbgs_iterator_poisson_matrix_gpu(mg_sdm[l].x, mg_a_poisson[l].coeff, mg_sdm[l].b,
                                 &mg_sdm[l], maxiteration, omega_sor, mg_sdm[l].is_aggregated);
                    // timer_stamp(15,STAMP_smooth);
                    // timer_stamp(6,STAMP_LEVEL);
                }

                // Prolongation to actual solution
                // timer_stamp0(STAMP_COMM_NEIGHBOR);
                geometry_halocell_update_selectively_gpu(mg_sdm[1].x, &mg_sdm[1], mg_sdm[1].is_aggregated);
                // timer_stamp(11,STAMP_COMM_NEIGHBOR);

                // timer_stamp0(STAMP_COMP);
                // timer_stamp0(STAMP_prolongation);
                multigrid_prolongation_linear_on_nonuniform_grid_gpu(rsd, mg_sdm[1].x, sdm, &mg_sdm[1], 0);
                // timer_stamp(17,STAMP_prolongation);
                // timer_stamp0(STAMP_LEVEL);
                nx = sdm->nx;
                ny = sdm->ny;
                nz = sdm->nz;
                    size_t size = (size_t)(nx + 2) * (size_t)(ny + 2) * (size_t)(nz + 2);
                adaptive_vector_add_gpu(sol, rsd, (size_t)size);
                // timer_stamp(10,STAMP_COMP);
                // timer_stamp0(STAMP_smooth);
                rbgs_iterator_poisson_matrix_gpu(sol, a_poisson->coeff, rhs, sdm, maxiteration, omega_sor, sdm->is_aggregated);
                // timer_stamp(15,STAMP_smooth);
                // timer_stamp(5,STAMP_LEVEL);
                
                // timer_stamp0(STAMP_residual);
                multigrid_residual_gpu(rsd, a_poisson->coeff, sol, rhs, sdm, sdm->is_aggregated);
                // timer_stamp(18,STAMP_residual);

                vv_dot_3d_matrix_gpu(&rsd_val, rsd, rsd, nx, ny, nz, sdm->is_aggregated);

                // if (myrank == 0) {
                //     printf("[MG] cycle = %d, Error: %e %e %e\n", cyc, sqrt(rsd_val), sqrt(res0tol), sqrt(rsd_val/res0tol));
                // }

                if (sqrt(rsd_val/res0tol) < tolerance) break;


            }

            // if(myrank==0) printf("[MG] Total %d V-cycles end\n", cyc);

            CUDA_CHECK(cudaFree(rsd));


            break;

        default:
            if (myrank == 0) printf("[Error] level_case should be 1, 2, or 3. Current: %d\n", level_case);
            MPI_Finalize();
            exit(EXIT_FAILURE);
    }



    


    if (lv_aggregation_x != 0) {
        MPI_Type_free(&ddtype_gatherv_x);
        MPI_Type_free(&ddtype_send_x);
    }
    if (lv_aggregation_y != 0) {
        MPI_Type_free(&ddtype_gatherv_y);
        MPI_Type_free(&ddtype_send_y);
    }
    if (lv_aggregation_z != 0) {
        MPI_Type_free(&ddtype_gatherv_z);
        MPI_Type_free(&ddtype_send_z);
    }
    free(disps_gatherv_z);
    free(cnt_gatherv_z);
    free(disps_gatherv_y);
    free(cnt_gatherv_y);
    free(disps_gatherv_x);
    free(cnt_gatherv_x);
    free(cart_coord);
}





void multigrid_destroy_gpu()
{
    multigrid_workspace_release_gpu();

    if (mg_sdm != NULL) {
        for (int l = 1; l <= n_levels; ++l) {
            geometry_subdomain_destroy_gpu(&mg_sdm[l]);
            geometry_subdomain_ddt_destroy(&mg_sdm[l]);
        }

        free(mg_sdm);
        mg_sdm = NULL;
    }

    if (mg_a_poisson != NULL) {
        for (int l = 1; l <= n_levels; ++l) {
            matrix_poisson_destroy_gpu(&mg_a_poisson[l]);
        }

        free(mg_a_poisson);
        mg_a_poisson = NULL;
    }
}



