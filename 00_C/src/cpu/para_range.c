#include <stdio.h>

void para_range(int n1, int n2, int nprocs, int myrank, int *ista, int *iend)
{
    int iwork1, iwork2;

    iwork1 = (n2 - n1 + 1) / nprocs;
    iwork2 = (n2 - n1 + 1) % nprocs;

    *ista = myrank * iwork1 + n1 + (myrank < iwork2 ? myrank : iwork2);
    *iend = *ista + iwork1 - 1;
    if (iwork2 > myrank) {
        *iend += 1;
    }
}
