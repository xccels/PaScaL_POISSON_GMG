#include "timer.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Fortran: double precision, private :: t_zero(*), t_curr */
static double t_zero[STAMP_MAX];
static double t_curr;

/* Fortran: double precision, public :: t_array(64), t_array_reduce(64) */
double t_array[TIMER_MAX];
double t_array_reduce[TIMER_MAX];

/* Fortran: integer(kind=4) :: ntimer */
int ntimer = 0;

/* Fortran: character(len=64) :: t_str(64) */
char t_str[TIMER_MAX][TIMER_STRLEN + 1];

static void die_with_mpi_finalize(const char *msg)
{
    fprintf(stderr, "%s\n", msg);
    /* C MPI: MPI_Finalize takes no args */
    MPI_Finalize();
    exit(EXIT_FAILURE);
}

void timer_init(int n, const char str[TIMER_MAX][TIMER_STRLEN + 1])
{
    if (n > TIMER_MAX) {
        die_with_mpi_finalize("[Error] Maximum number of timer is 64");
    }
    if (n < 0) {
        die_with_mpi_finalize("[Error] timer_init: n must be >= 0");
    }

    ntimer = n;

    /* t_array(:)=0, t_array_reduce(:)=0, t_str(:)='null' */
    for (int i = 0; i < TIMER_MAX; ++i) {
        t_array[i] = 0.0;
        t_array_reduce[i] = 0.0;
        strncpy(t_str[i], "null", TIMER_STRLEN);
        t_str[i][TIMER_STRLEN] = '\0';
    }

    /* copy the first ntimer strings */
    for (int i = 0; i < ntimer; ++i) {
        /* ensure null-termination */
        strncpy(t_str[i], str[i], TIMER_STRLEN);
        t_str[i][TIMER_STRLEN] = '\0';
    }

    /* init stamps */
    for (int s = 0; s < STAMP_MAX; ++s) {
        t_zero[s] = 0.0;
    }
    t_curr = 0.0;
}

void timer_stamp0(int stamp_id)
{
    if (stamp_id < 0 || stamp_id >= STAMP_MAX) {
        die_with_mpi_finalize("[Error] timer_stamp0: stamp_id out of range");
    }
    t_zero[stamp_id] = MPI_Wtime();
}

void timer_stamp(int timer_id, int stamp_id)
{
    /* Fortran timer_id is 1-based; in C we keep it 1-based externally too,
       so store in index timer_id-1. */
    if (timer_id < 1 || timer_id > TIMER_MAX) {
        die_with_mpi_finalize("[Error] timer_stamp: timer_id out of range");
    }
    if (stamp_id < 0 || stamp_id >= STAMP_MAX) {
        die_with_mpi_finalize("[Error] timer_stamp: stamp_id out of range");
    }

    const int idx = timer_id - 1;

    t_curr = MPI_Wtime();
    t_array[idx] += (t_curr - t_zero[stamp_id]);
    t_zero[stamp_id] = t_curr;
}

void timer_start(int timer_id)
{
    if (timer_id < 1 || timer_id > TIMER_MAX) {
        die_with_mpi_finalize("[Error] timer_start: timer_id out of range");
    }
    t_array[timer_id - 1] = MPI_Wtime();
}

void timer_end(int timer_id)
{
    if (timer_id < 1 || timer_id > TIMER_MAX) {
        die_with_mpi_finalize("[Error] timer_end: timer_id out of range");
    }
    const int idx = timer_id - 1;
    t_array[idx] = MPI_Wtime() - t_array[idx];
}

double timer_elapsed(int timer_id)
{
    if (timer_id < 1 || timer_id > TIMER_MAX) {
        die_with_mpi_finalize("[Error] timer_elapsed: timer_id out of range");
    }
    return MPI_Wtime() - t_array[timer_id - 1];
}

void timer_reduction(MPI_Comm comm)
{
    for (int i = 0; i < TIMER_MAX; ++i) {
        t_array_reduce[i] = 0.0;
    }

    /* C MPI: returns int error code */
    int ierr = MPI_Reduce(t_array, t_array_reduce, ntimer,
                          MPI_DOUBLE, MPI_SUM, 0, comm);

    if (ierr != MPI_SUCCESS) {
        die_with_mpi_finalize("[Error] MPI_Reduce failed in timer_reduction");
    }
}

void timer_output(int myrank, int nprocs)
{
    if (myrank == 0) {
        for (int i = 0; i < ntimer; ++i) {
            /* mimic Fortran: if(trim(t_str(i)).ne.'null') */
            if (strncmp(t_str[i], "null", 4) != 0) {
                /* print average across ranks: t_array_reduce(i)/nprocs */
                printf("[Timer] %-35s : (%3d) : %16.9f\n",
                       t_str[i], i + 1, t_array_reduce[i] / (double)nprocs);
            }
        }
    }
}
