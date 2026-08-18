#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime.h>
// #include "timer.h"
#include "geometry.h"
#include "matrix.h"
#include "multigrid.h"
#include "mpi_topology.h"
#include "para_range.h"
#include "rbgs_poisson_matrix.h"
// #include "global.h"

#define PI 3.14159265358979323846
#define IDX(i,j,k,ni,nj) ((i)*(ni) + (j)*(nj) + (k))

#define CUDA_CHECK(call) do {                                      \
    cudaError_t err = (call);                                      \
    if (err != cudaSuccess) {                                      \
        fprintf(stderr, "CUDA error at %s:%d: %s\n",               \
                __FILE__, __LINE__, cudaGetErrorString(err));      \
        exit(EXIT_FAILURE);                                        \
    }                                                             \
} while (0)

extern "C" {
// int nprocs = 1;
// int myrank = 0;
int lv_gdm_coarsest_x = 999;
int lv_gdm_coarsest_y = 999;
int lv_gdm_coarsest_z = 999;
}



__global__
void initialize_rhs_and_reference_kernel(double *d_x, double *d_b, double *d_ref,
                                         const double *xg, const double *yg, const double *zg,
                                         int nx, int ny, int nz, int ni, int nj)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i > nx + 1 || j > ny + 1 || k > nz + 1) return;

    int idx = IDX(i, j, k, ni, nj);

    double ref = cos(xg[i] * PI) * cos(yg[j] * PI) * cos(zg[k] * PI);

    d_x[idx]   = 0.0;
    d_ref[idx] = ref;

    if (i >= 1 && i <= nx &&
        j >= 1 && j <= ny &&
        k >= 1 && k <= nz)
    {
        d_b[idx] = -3.0 * PI * PI * ref;
    }
    else
    {
        d_b[idx] = 0.0;
    }
}


__global__
void compute_rms_kernel(const double *d_x, const double *d_ref, double *d_sum,
                        int nx, int ny, int nz, int ni, int nj)
{
    extern __shared__ double s_sum[];

    int local_tid = threadIdx.x;
    int global_tid = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    int ncell = nx * ny * nz;

    double local_sum = 0.0;

    // 每个线程处理一个或多个网格点
    for (int p = global_tid; p < ncell; p += stride) {
        int k = p % nz + 1;
        int j = (p / nz) % ny + 1;
        int i = p / (ny * nz) + 1;

        int idx = IDX(i, j, k, ni, nj);
        double error = d_x[idx] - d_ref[idx];

        local_sum += error * error;
    }

    // block 内归约
    s_sum[local_tid] = local_sum;
    __syncthreads();

    for (int offset = blockDim.x / 2; offset > 0; offset >>= 1) {
        if (local_tid < offset) {
            s_sum[local_tid] += s_sum[local_tid + offset];
        }
        __syncthreads();
    }

    // 每个 block 只执行一次原子加法
    if (local_tid == 0) {
        atomicAdd(d_sum, s_sum[0]);
    }
}

static void read_input_file(const char *filename,
                            int *nx, int *ny, int *nz,
                            double *ox, double *oy, double *oz,
                            double *lx, double *ly, double *lz,
                            double *ax, double *ay, double *az,
                            int *npx, int *npy, int *npz,
                            int *maxiteration, double *tolerance,
                            int *number_of_vcycles, int *number_of_levels,
                            int *aggregation_method, int *aggregation_level,
                            double *omega_sor)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Cannot open input file: %s\n", filename);
        exit(EXIT_FAILURE);
    }
    char line[256], key[128];
    int ival;
    double dval;
    while (fgets(line, sizeof(line), fp)) {
        if (sscanf(line, "%127s %lf", key, &dval) == 2) {
            if      (strcmp(key, "ox") == 0) *ox = dval;
            else if (strcmp(key, "oy") == 0) *oy = dval;
            else if (strcmp(key, "oz") == 0) *oz = dval;
            else if (strcmp(key, "lx") == 0) *lx = dval;
            else if (strcmp(key, "ly") == 0) *ly = dval;
            else if (strcmp(key, "lz") == 0) *lz = dval;
            else if (strcmp(key, "ax") == 0) *ax = dval;
            else if (strcmp(key, "ay") == 0) *ay = dval;
            else if (strcmp(key, "az") == 0) *az = dval;
            else if (strcmp(key, "tolerance") == 0) *tolerance = dval;
            else if (strcmp(key, "omega_sor") == 0) *omega_sor = dval;
        }
        if (sscanf(line, "%127s %d", key, &ival) == 2) {
            if      (strcmp(key, "nx") == 0) *nx = ival;
            else if (strcmp(key, "ny") == 0) *ny = ival;
            else if (strcmp(key, "nz") == 0) *nz = ival;
            else if (strcmp(key, "npx") == 0) *npx = ival;
            else if (strcmp(key, "npy") == 0) *npy = ival;
            else if (strcmp(key, "npz") == 0) *npz = ival;
            else if (strcmp(key, "maxiteration") == 0) *maxiteration = ival;
            else if (strcmp(key, "number_of_vcycles") == 0) *number_of_vcycles = ival;
            else if (strcmp(key, "number_of_levels") == 0) *number_of_levels = ival;
            else if (strcmp(key, "aggregation_method") == 0) *aggregation_method = ival;
            else if (strcmp(key, "aggregation_level") == 0) *aggregation_level = ival;
        }
    }
    fclose(fp);
    // *aggregation_method = 0;
    // *aggregation_level = 0;
}

static void output_solution_gpu(const subdomain *sdm, const double *d_ref, int myrank)
{
    int nx = sdm->nx;
    int ny = sdm->ny;
    int nz = sdm->nz;

    int ni = (ny + 2) * (nz + 2);
    int nj = nz + 2;

    size_t size3d =
        (size_t)(nx + 2) *
        (size_t)(ny + 2) *
        (size_t)(nz + 2);

    double *h_x   = (double*)malloc(size3d * sizeof(double));
    double *h_b   = (double*)malloc(size3d * sizeof(double));
    double *h_ref = (double*)malloc(size3d * sizeof(double));

    double *h_xg = (double*)malloc((nx + 2) * sizeof(double));
    double *h_yg = (double*)malloc((ny + 2) * sizeof(double));
    double *h_zg = (double*)malloc((nz + 2) * sizeof(double));

    CUDA_CHECK(cudaMemcpy(h_x,   sdm->x, size3d * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_b,   sdm->b, size3d * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_ref, d_ref,  size3d * sizeof(double), cudaMemcpyDeviceToHost));

    CUDA_CHECK(cudaMemcpy(h_xg, sdm->xg, (nx + 2) * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_yg, sdm->yg, (ny + 2) * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_zg, sdm->zg, (nz + 2) * sizeof(double), cudaMemcpyDeviceToHost));

    char filename[256];
    // snprintf(filename, sizeof(filename), "./result/solution.%03d", myrank);
    snprintf(filename, sizeof(filename), "./result/solution.plt");

    FILE *fp = fopen(filename, "w");
    if (!fp) {
        printf("Cannot open output file: %s\n", filename);
        exit(1);
    }

    for (int k = 0; k <= nz + 1; k++) {
        for (int j = 0; j <= ny + 1; j++) {
            for (int i = 0; i <= nx + 1; i++) {
                int idx = IDX(i, j, k, ni, nj);

                fprintf(fp, "% .16e % .16e % .16e % .16e % .16e % .16e\n",
                        h_xg[i],
                        h_yg[j],
                        h_zg[k],
                        h_x[idx],
                        h_ref[idx],
                        h_b[idx]);
            }
        }
    }

    fclose(fp);

    free(h_x);
    free(h_b);
    free(h_ref);

    free(h_xg);
    free(h_yg);
    free(h_zg);
}


int main(int argc, char **argv)
{
    int npx, npy, npz;
    int nx, ny, nz;
    double ox, oy, oz;
    double lx, ly, lz;
    double ax, ay, az;
    double alpha_x, alpha_y, alpha_z;
    // aggretation method 
    // 0 : no aggregation
    // 1 : single aggregation
    // 2 : adaptive aggregation (not implemented yet)
    int i, j, k, idx, maxiteration, number_of_vcycles, number_of_levels, aggregation_method, aggregation_level;
    double tolerance, omega_sor;


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int ngpus = 0;
    cudaGetDeviceCount(&ngpus);
    cudaSetDevice(myrank % ngpus);



    
    // 读取文件 PARA_INPUT.inp
    // namelist /meshes/ nx, ny, nz
    // namelist /origin/ ox, oy, oz
    // namelist /length/ lx, ly, lz
    // namelist /mesh_stretch/ ax, ay, az
    // namelist /procs/ npx, npy, npz
    // namelist /control/ maxiteration, tolerance, number_of_vcycles, number_of_levels, aggregation_method, aggregation_level, omega_sor
    // namelist /coefficients/ alpha_x, alpha_y, alpha_z

    const char *input_file = (argc >= 2) ? argv[1] : "run/PARA_INPUT.inp";

    read_input_file(input_file, &nx, &ny, &nz, &ox, &oy, &oz, &lx, &ly, &lz,
                    &ax, &ay, &az, &npx, &npy, &npz,
                    &maxiteration, &tolerance,
                    &number_of_vcycles, &number_of_levels,
                    &aggregation_method, &aggregation_level, &omega_sor);

    if(myrank==0)
    {
        printf("[GPU GMG] nx=%d ny=%d nz=%d\n", nx, ny, nz);
        printf("[GPU GMG] levels=%d vcycles=%d maxiter=%d tol=%e omega=%f\n",
                number_of_levels, number_of_vcycles, maxiteration, tolerance, omega_sor);
    }


    // char timer_str[TIMER_MAX][TIMER_STRLEN + 1];
    // for (int i = 0; i < TIMER_MAX; ++i) {
    //     strncpy(timer_str[i], "null", TIMER_STRLEN);
    //     timer_str[i][TIMER_STRLEN] = '\0';
    // }

    // strcpy(timer_str[0],  "[Main] multigrid_create             ");
    // strcpy(timer_str[1],  "[Main] multigrid_solve_vcycle       ");
    // strcpy(timer_str[2],  "[Main]            ");
    // strcpy(timer_str[3],  "[Main]                ");
    // strcpy(timer_str[4],  "[level] finest                      ");
    // strcpy(timer_str[5],  "[level] intermediate                ");
    // strcpy(timer_str[6],  "[level] coarest                     ");
    // strcpy(timer_str[7],  "[level] aggregation                 ");
    // strcpy(timer_str[8],  "[level]                     ");
    // strcpy(timer_str[9],  "[com] computation                   ");
    // strcpy(timer_str[10], "[com] comm_neighbor                 ");
    // strcpy(timer_str[11], "[com] comm_Allreduce                ");
    // strcpy(timer_str[12], "[com] Gather                        ");
    // strcpy(timer_str[13], "[com] Scatter                       ");

    // strcpy(timer_str[14], "[divide] smooth                     ");
    // strcpy(timer_str[15], "[divide] restriction                ");
    // strcpy(timer_str[16], "[divide] prolongation               ");
    // strcpy(timer_str[17], "[divide] residual                   ");









    int times;



    for(times=1; times<=20; times++)
    {


    // timer_init(18, timer_str);
    // cudaEvent_t ev_start, ev_stop;
    // CUDA_CHECK(cudaEventCreate(&ev_start));
    // CUDA_CHECK(cudaEventCreate(&ev_stop));
    // CUDA_CHECK(cudaEventRecord(ev_start));


    np_dim[0] = npx; np_dim[1] = npy; np_dim[2] = npz;
    period[0] = period[1] = period[2] = 0;
    mpi_topology_create();
    mpi_boundary_create();


    domain g_domain;
    subdomain s_domain;
    matrix_poisson a_poisson;
    memset(&g_domain, 0, sizeof(domain));
    memset(&s_domain, 0, sizeof(subdomain));
    memset(&a_poisson, 0, sizeof(matrix_poisson));

    

    geometry_domain_create_gpu(&g_domain, nx, ny, nz, ox, oy, oz, lx, ly, lz, ax, ay, az, period);
    // if (myrank == 0) {
    //         printf("[Poisson] Geometry and matrix size initialized.\n");
    //     }
    geometry_subdomain_create_gpu(&s_domain, &g_domain);
    geometry_subdomain_ddt_create(&s_domain);

    int ni = (s_domain.ny + 2) * (s_domain.nz + 2);
    int nj = (s_domain.nz + 2);
    size_t size3d = (size_t)(s_domain.nx + 2) * (s_domain.ny + 2) * (s_domain.nz + 2);

    double *d_ref = NULL;
    CUDA_CHECK(cudaMalloc((void**)&d_ref, size3d * sizeof(double)));
    CUDA_CHECK(cudaMemset(d_ref, 0, size3d * sizeof(double)));

    dim3 block3d(8, 8, 4);
    dim3 grid3d((s_domain.nx + 2 +block3d.x - 1) / block3d.x,
                (s_domain.ny + 2 +block3d.y - 1) / block3d.y,
                (s_domain.nz + 2 +block3d.z - 1) / block3d.z);
    initialize_rhs_and_reference_kernel<<<grid3d, block3d>>>(
        s_domain.x, s_domain.b, d_ref, s_domain.xg, s_domain.yg, s_domain.zg,
        s_domain.nx, s_domain.ny, s_domain.nz, ni, nj);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

        // if (myrank == 0) {
        //     printf("[Poisson] Geometry and rhs constructed.\n");
        // }




    // 主体部分

    MPI_Barrier(MPI_COMM_WORLD);
    double setup_start = MPI_Wtime();

    
    // timer_stamp0(STAMP_COMP);
    matrix_poisson_create_gpu(&a_poisson, &s_domain);
    // timer_stamp(10,STAMP_COMP);

    // if (myrank == 0) {
    //     printf("[Poisson] Poisson matrix constructed.\n");
    // }
    // if (myrank == 0) {
    //     printf("[Poisson] Start solving equations.\n");
    // }
    
    // double t0 = MPI_Wtime();

    // timer_stamp0(STAMP_MAIN);
    multigrid_create_gpu(&s_domain, number_of_levels, number_of_vcycles, aggregation_method, aggregation_level);
    // timer_stamp(1, STAMP_MAIN);

    CUDA_CHECK(cudaDeviceSynchronize());
    double setup_end = MPI_Wtime();
    double setup_local = setup_end - setup_start;





    MPI_Barrier(MPI_COMM_WORLD);
    double solve_start = MPI_Wtime();

    multigrid_solve_vcycle_gpu(s_domain.x, &a_poisson, s_domain.b, &s_domain, maxiteration, tolerance, omega_sor);
    // timer_stamp(2, STAMP_MAIN);
    CUDA_CHECK(cudaDeviceSynchronize());
    double solve_end = MPI_Wtime();
    double solve_local = solve_end - solve_start;

    multigrid_destroy_gpu();

   
    double setup_global = 0.0;
    double solve_global = 0.0;
    MPI_Reduce(&setup_local, &setup_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&solve_local, &solve_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

// if (myrank == 0) {
//     printf("Setup time = %.6f s\n", setup_global);
//     printf("Solve time = %.6f s\n", solve_global);
// }
    

   

    double *d_rms_sum = NULL;
    double h_rms_sum = 0.0;

    CUDA_CHECK(cudaMalloc((void**)&d_rms_sum, sizeof(double)));

    CUDA_CHECK(cudaMemset(d_rms_sum, 0, sizeof(double)));

    int ncell = s_domain.nx * s_domain.ny * s_domain.nz;
    int block1d = 256;
    int grid1d = (ncell + block1d - 1) / block1d;

    if (grid1d > 4096) grid1d = 4096;

    size_t shared_bytes = block1d * sizeof(double);

    // timer_stamp0(STAMP_COMP);
    compute_rms_kernel<<<grid1d, block1d, shared_bytes>>>(s_domain.x, d_ref, d_rms_sum, s_domain.nx, s_domain.ny, s_domain.nz, ni, nj);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    // timer_stamp(10,STAMP_COMP);

    CUDA_CHECK(cudaMemcpy(&h_rms_sum, d_rms_sum, sizeof(double), cudaMemcpyDeviceToHost));

    // CUDA_CHECK(cudaEventRecord(ev_stop));
    // CUDA_CHECK(cudaEventSynchronize(ev_stop));
    // float elapsed_ms = 0.0f;
    // CUDA_CHECK(cudaEventElapsedTime(&elapsed_ms, ev_start, ev_stop));


    double h_rms_global = 0.0;

    // timer_stamp0(STAMP_COMM_ALLREDUCE);
    MPI_Allreduce(&h_rms_sum, &h_rms_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    // timer_stamp(12,STAMP_COMM_ALLREDUCE);

    if (myrank == 0) {
        printf("[GPU GMG] RMS = %e\n", h_rms_global / (double)(nx * ny * nz));
        printf("[GPU GMG] Setup time = %.6f s\n", setup_global);
        printf("[GPU GMG] Solve time = %.6f s\n", solve_global);
    }

    // timer_reduction(MPI_COMM_WORLD);
    // timer_output(myrank, nprocs);
   

    if (myrank == 0) {
        system("mkdir -p ./result");
    }
    // output_solution_gpu(&s_domain, d_ref, myrank);

    CUDA_CHECK(cudaFree(d_rms_sum));
    CUDA_CHECK(cudaFree(d_ref));

    matrix_poisson_destroy_gpu(&a_poisson);

    mpi_topology_destroy();
    geometry_subdomain_ddt_destroy(&s_domain);
    geometry_subdomain_destroy_gpu(&s_domain);
    geometry_domain_destroy_gpu(&g_domain);
    // CUDA_CHECK(cudaEventDestroy(ev_start));
    // CUDA_CHECK(cudaEventDestroy(ev_stop));

    if (myrank == 0) {
        printf("[GPU GMG] Memory deallocated.\n");
    }




    }











    MPI_Finalize();
    return 0;
}
