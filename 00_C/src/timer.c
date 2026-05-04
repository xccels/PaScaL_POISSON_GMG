// timer.c
#include "timer.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#define MAX_TIMERS 64

static double t_zero[MAX_TIMERS], t_curr;
double t_array[MAX_TIMERS], t_array_reduce[MAX_TIMERS];
int ntimer;
char t_str[MAX_TIMERS][64];

void timer_init(int n, char str[][64]) {
    int i;
    if (n > MAX_TIMERS) {
        fprintf(stderr, "[Error] Maximum number of timers is %d\n", MAX_TIMERS);
        int ierr;
        MPI_Finalize();
        exit(1);
    }

    ntimer = n;
    for (i = 0; i < ntimer; ++i) {
        t_array[i] = 0.0;
        t_array_reduce[i] = 0.0;
        memset(t_str[i], 0, 64);               // 🔧 add this
        strncpy(t_str[i], str[i], 63);         // ✅ secure copy, leave space for \0
        t_str[i][63] = '\0';  
    }
}

void timer_stamp0(int stamp_id) {
    t_zero[stamp_id] = MPI_Wtime();
}

void timer_stamp(int timer_id, int stamp_id) {
    t_curr = MPI_Wtime();
    t_array[timer_id] += t_curr - t_zero[stamp_id];
    t_zero[stamp_id] = t_curr;
}

void timer_start(int timer_id) {
    t_array[timer_id] = MPI_Wtime();
}

void timer_end(int timer_id) {
    t_array[timer_id] = MPI_Wtime() - t_array[timer_id];
}

double timer_elapsed(int timer_id) {
    return MPI_Wtime() - t_array[timer_id];
}

void timer_reduction() {
    int ierr;
    for (int i = 0; i < ntimer; ++i) t_array_reduce[i] = 0.0;
    MPI_Reduce(t_array, t_array_reduce, ntimer, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}

void timer_output(int myrank, int nprocs) {
    if (myrank == 0) {
        for (int i = 0; i < ntimer; ++i) {
            if (strcmp(t_str[i], "null") != 0) {
                printf("[Timer] %-35s : (%2d) : %16.9f\n", t_str[i], i + 1, t_array_reduce[i] / nprocs);
            }
        }
    }
}
