#ifndef TIMER_H
#define TIMER_H

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Match Fortran: t_array(64), t_array_reduce(64), t_str(64), ntimer */
#define TIMER_MAX 64
#define TIMER_STRLEN 64

/* Your stamp ids in main go up to 14, so stamp array must be >= 14.
   For safety, keep it 64 (same as TIMER_MAX). */
#define STAMP_MAX 64

/* Public state (like Fortran module public variables) */
extern double t_array[TIMER_MAX];
extern double t_array_reduce[TIMER_MAX];
extern int    ntimer;
extern char   t_str[TIMER_MAX][TIMER_STRLEN + 1];

/* API (like Fortran contained subroutines/functions) */
void timer_init(int n, const char str[TIMER_MAX][TIMER_STRLEN + 1]);
void timer_stamp0(int stamp_id);
void timer_stamp(int timer_id, int stamp_id);
void timer_start(int timer_id);
void timer_end(int timer_id);
double timer_elapsed(int timer_id);
void timer_reduction(MPI_Comm comm);
void timer_output(int myrank, int nprocs);

/* Optional: define your stamp IDs in one place (same as Fortran parameters) */
enum {
    STAMP_MAIN            = 1,
    STAMP_LEVEL           = 5,
    STAMP_AGG             = 8,
    STAMP_COMP            = 10,
    STAMP_COMM_NEIGHBOR   = 11,
    STAMP_COMM_ALLREDUCE  = 12,
    STAMP_COMM_GATHER     = 13,
    STAMP_COMM_SCATTER    = 14,

    STAMP_smooth          = 15,
    STAMP_restriction     = 16,
    STAMP_prolongation    = 17,
    STAMP_residual        = 18
};

#ifdef __cplusplus
}
#endif

#endif /* TIMER_H */
