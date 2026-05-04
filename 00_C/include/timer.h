#ifndef TIMER_H
#define TIMER_H

#define MAX_TIMERS 64

extern double t_array[MAX_TIMERS];
extern double t_array_reduce[MAX_TIMERS];
extern int ntimer;
extern char t_str[MAX_TIMERS][64];

void timer_init(int n, char str[][64]);
void timer_stamp0(int stamp_id);
void timer_stamp(int timer_id, int stamp_id);
void timer_start(int timer_id);
void timer_end(int timer_id);
double timer_elapsed(int timer_id);
void timer_reduction();
void timer_output(int myrank, int nprocs);

#endif
