#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <mpi.h>
#include <stdbool.h>
#include <mpi_topology.h>
#include <para_range.h>
#include "global.h"

typedef struct {
    int nx, ny, nz;
    double lx, ly, lz;
    double ox, oy, oz;
    double *dxm, *dym, *dzm;
    double *dxg, *dyg, *dzg;
    double *xg, *yg, *zg;
    int is_periodic[3];
} domain;

typedef struct {
    int nx, ny, nz;
    double lx, ly, lz;
    double ox, oy, oz;
    double *dxm, *dym, *dzm;
    double *dxg, *dyg, *dzg;
    double *xg, *yg, *zg;
    double *x, *b, *r;

    int is_periodic[3];
    int is_aggregated[3];

    int ista, iend, jsta, jend, ksta, kend;
    // int ddt_yz_plane_x0, ddt_yz_plane_x1, ddt_yz_plane_xn, ddt_yz_plane_xn1;
    // int ddt_xz_plane_y0, ddt_xz_plane_y1, ddt_xz_plane_yn, ddt_xz_plane_yn1;
    // int ddt_xy_plane_z0, ddt_xy_plane_z1, ddt_xy_plane_zn, ddt_xy_plane_zn1;
    // int ddt_inner_domain;
    MPI_Datatype ddt_yz_plane_x0, ddt_yz_plane_x1, ddt_yz_plane_xn, ddt_yz_plane_xn1;
    MPI_Datatype ddt_xz_plane_y0, ddt_xz_plane_y1, ddt_xz_plane_yn, ddt_xz_plane_yn1;
    MPI_Datatype ddt_xy_plane_z0, ddt_xy_plane_z1, ddt_xy_plane_zn, ddt_xy_plane_zn1;
    MPI_Datatype ddt_inner_domain;
    bool is_x0_boundary, is_x1_boundary;
    bool is_y0_boundary, is_y1_boundary;
    bool is_z0_boundary, is_z1_boundary;
} subdomain;

void geometry_domain_create(domain* dom, int nx, int ny, int nz,
                            double lx, double ly, double lz,
                            double ox, double oy, double oz,
                            double ax, double ay, double az,
                            int is_periodic[3]);
void geometry_domain_destroy(domain* dom);

void geometry_subdomain_create(subdomain* sdm, const domain* gdm);
void geometry_subdomain_destroy(subdomain* sdm);

void geometry_subdomain_ddt_create(subdomain *sdm);
void geometry_subdomain_ddt_destroy(subdomain *sdm);

void geometry_halocell_update(double ***u, subdomain *sdm);
void geometry_halocell_update_selectively(double *u, subdomain *sdm, int *is_serial);

void geometry_boundary_value_calculate(double ***u, subdomain *sdm, domain *gdm, boundary_comm *comm_boundary);
void geometry_boundary_normal_derivative(double ***dudx, double ***dudy, double ***dudz, double ***u, subdomain *sdm);
void geometry_boundary_surface_integral(double ***u, double ***dudx, double ***dudy, double ***dudz,subdomain *sdm, domain *gdm);

#endif