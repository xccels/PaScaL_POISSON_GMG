#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "geometry.h"     
#include "matrix.h"
#include "multigrid.h"
#include "mpi_topology.h"
#include "para_range.h"
#include "rbgs_poisson_matrix.h" 
#include "global.h"    

#define PI 3.14159265358979323846

#define IDX(i,j,k, nx,ny,nz) ((i)*(ny+2)*(nz+2) + (j)*(nz+2) + (k))

int nprocs, myrank;

void output(double *phi, double *xg, double *yg)
{
   int i, j, idx;

   // Write results for the whole domain 
   FILE *file = fopen("phi.plt", "w");
   fprintf(file, "VARIABLES=X, Y, phi \n");
   fprintf(file, "zone t=\"%d\" ",1);
   fprintf(file, "i=%d ",256);
   fprintf(file, "j=%d\n",256);

   for (i = 1; i <= 256; i++)
   {
      for (j = 1; j <= 256; j++)
      {
        idx = IDX(i,j,128,256,256,256);
        fprintf(file, "%.6f  ", xg[i]);
        fprintf(file, "%.6f  ", yg[j]);
        fprintf(file, "%e  ", phi[idx]);
        fprintf(file, "\n");
      }
   }

   fclose(file);
}

int main(int argc, char** argv) 
{
    

    int npx, npy, npz;
    int nx, ny, nz;
    double ox, oy, oz;
    double lx, ly, lz;
    double ax, ay, az;
    double alpha_x, alpha_y, alpha_z;

    domain g_domain;
    subdomain s_domain;
    matrix_poisson a_poisson;

    double *ref_sub; 
    double u_x, u_y, u_z;

    double rms, rms_local;

    // aggretation method 
    // 0 : no aggregation
    // 1 : single aggregation
    // 2 : adaptive aggregation (not implemented yet)
    int i, j, k, idx, maxiteration, number_of_vcycles, number_of_levels, aggregation_method, aggregation_level;
    double tolerance, t0, omega_sor;

    int times;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // 读取文件 PARA_INPUT.inp
    // namelist /meshes/ nx, ny, nz
    // namelist /origin/ ox, oy, oz
    // namelist /length/ lx, ly, lz
    // namelist /mesh_stretch/ ax, ay, az
    // namelist /procs/ npx, npy, npz
    // namelist /control/ maxiteration, tolerance, number_of_vcycles, number_of_levels, aggregation_method, aggregation_level, omega_sor
    // namelist /coefficients/ alpha_x, alpha_y, alpha_z
    FILE *fp = fopen("run/PARA_INPUT.inp", "r");
    if (fp == NULL) {
        printf("Cannot open file PARA_INPUT.inp\n");
        return 1; 
    }

    char line[40];
    char key[40];
    int ival;
    double dval;
    while (fgets(line, sizeof(line), fp) != NULL) 
    {
        if (sscanf(line, "%s %lf", key, &dval) == 2)
        {
            if (strcmp(key, "ox") == 0) ox = dval;
            else if (strcmp(key, "oy") == 0) oy = dval;
            else if (strcmp(key, "oz") == 0) oz = dval;
            else if (strcmp(key, "lx") == 0) lx = dval;
            else if (strcmp(key, "ly") == 0) ly = dval;
            else if (strcmp(key, "lz") == 0) lz = dval;
            else if (strcmp(key, "ax") == 0) ax = dval;
            else if (strcmp(key, "ay") == 0) ay = dval;
            else if (strcmp(key, "az") == 0) az = dval;
            else if (strcmp(key, "alpha_x") == 0) alpha_x = dval;
            else if (strcmp(key, "alpha_y") == 0) alpha_y = dval;
            else if (strcmp(key, "alpha_z") == 0) alpha_z = dval;
            else if (strcmp(key, "tolerance") == 0) tolerance = dval;
            else if (strcmp(key, "omega_sor") == 0) omega_sor = dval;
        }
        if (sscanf(line, "%s %d", key, &ival) == 2) 
        {
            if (strcmp(key, "nx") == 0) nx = ival;
            else if (strcmp(key, "ny") == 0) ny = ival;
            else if (strcmp(key, "nz") == 0) nz = ival;
            else if (strcmp(key, "npx") == 0) npx = ival;
            else if (strcmp(key, "npy") == 0) npy = ival;
            else if (strcmp(key, "npz") == 0) npz = ival;
            else if (strcmp(key, "maxiteration") == 0) maxiteration = ival;
            else if (strcmp(key, "number_of_vcycles") == 0) number_of_vcycles = ival;
            else if (strcmp(key, "number_of_levels") == 0) number_of_levels = ival;
            else if (strcmp(key, "aggregation_method") == 0) aggregation_method = ival;
            else if (strcmp(key, "aggregation_level") == 0) aggregation_level = ival;
        }
    }
    fclose(fp);
    // printf("nx=%d ny=%d nz=%d\n", nx, ny, nz);
    // printf("ox=%lf oy=%lf oz=%lf\n", ox, oy, oz);
    // printf("lx=%lf ly=%lf lz=%lf\n", lx, ly, lz);
    // printf("ax=%lf ay=%lf az=%lf\n", ax, ay, az);
    // printf("npx=%d npy=%d npz=%d\n", npx, npy, npz);
    // printf("maxiteration=%d tolerance=%e\n", maxiteration, tolerance);
    // printf("number_of_vcycles=%d number_of_levels=%d\n", number_of_vcycles, number_of_levels);
    // printf("aggregation_method=%d aggregation_level=%d omega_sor=%lf\n", aggregation_method, aggregation_level, omega_sor);
    // printf("alpha_x=%lf alpha_y=%lf alpha_z=%lf\n", alpha_x, alpha_y, alpha_z);


    for(times=1; times<=1; times++)
    {
        np_dim[0] = npx; np_dim[1] = npy; np_dim[2] = npz;
        period[0] = period[1] = period[2] = 0;
        mpi_topology_create();
        mpi_boundary_create();
        geometry_domain_create(&g_domain, nx, ny, nz, ox, oy, oz, lx, ly, lz, ax, ay, az, period);
        if (myrank == 0) {
            printf("[Poisson] Geometry and matrix size initialized.\n");
        }
        geometry_subdomain_create(&s_domain, &g_domain);
        geometry_subdomain_ddt_create(&s_domain);

        size_t size3d = (s_domain.nx + 2) * (s_domain.ny + 2) * (s_domain.nz + 2);
        ref_sub = calloc(size3d, sizeof(double));

        for (i = 1; i <= s_domain.nx; i++)
        {
           for (j = 1; j <= s_domain.ny; j++)
           {
                for (k = 1; k <= s_domain.nz; k++)
                {
                    idx = IDX(i,j,k,s_domain.nx,s_domain.ny,s_domain.nz);
                    s_domain.x[idx] = 0.0;
                    s_domain.b[idx] = -cos(s_domain.xg[i] * PI) * cos(s_domain.yg[j] * PI) * cos(s_domain.zg[k] * PI) * 3.0 * PI * PI;
                    ref_sub[idx] = cos(s_domain.xg[i] * PI) * cos(s_domain.yg[j] * PI) * cos(s_domain.zg[k] * PI);
                }
           }
        }
        
        if (myrank == 0) {
            printf("[Poisson] Geometry and rhs constructed.\n");
        }

        matrix_poisson_create(&a_poisson, &s_domain);
        if (myrank == 0) {
            printf("[Poisson] Poisson matrix constructed.\n");
        }

        if (myrank == 0) {
            printf("[Poisson] Start solving equations.\n");
        }

        t0 = MPI_Wtime();

        multigrid_create(&s_domain, number_of_levels, number_of_vcycles, aggregation_method, aggregation_level);
        multigrid_solve_vcycle(s_domain.x, &a_poisson, s_domain.b, &s_domain, maxiteration, tolerance, omega_sor);
        multigrid_destroy();
        

        rms = 0.0;
        rms_local = 0.0;
        for (i = 1; i <= s_domain.nx; i++)
        {
           for (j = 1; j <= s_domain.ny; j++)
           {
                for (k = 1; k <= s_domain.nz; k++)
                {
                    // printf("[Error3] \n");
                    idx = IDX(i,j,k,s_domain.nx,s_domain.ny,s_domain.nz);
                    rms_local = rms_local + (s_domain.x[idx] - ref_sub[idx]) * (s_domain.x[idx] - ref_sub[idx]);
                }
           }
        }
        // printf("[Error1] \n");
        // output(s_domain.x, s_domain.xg, s_domain.yg);
        // printf("[Error2] \n");

        MPI_Allreduce(&rms_local, &rms, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        // printf("[Error3] \n");
        if (myrank == 0) {
            printf("[Poisson] RMS = %e\n", rms / (g_domain.nx * g_domain.ny * g_domain.nz));
            printf("[Poisson] Solution obtained. Execution time = %f\n", MPI_Wtime() - t0);
        }

        free(ref_sub);

        matrix_poisson_destroy(&a_poisson);
        
        mpi_topology_destroy();
        geometry_subdomain_ddt_destroy(&s_domain);
        geometry_subdomain_destroy(&s_domain);
        geometry_domain_destroy(&g_domain);

        if (myrank == 0) {
            printf("[Poisson] Memory deallocated.\n");
        }
        


    }
    

   
    MPI_Finalize();
    return 0;
}