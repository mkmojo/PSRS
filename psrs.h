#include <err.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "mpi.h"

#define VALIDATION_THRESHOLD 200000 /* if nInts < threshold, validate sort */

#ifdef DEBUG
#define DBPRINT(str) do { printf str ; fflush(stdout);} while(0)
#else
#define DBPRINT(str) do {} while (0)
#endif

#define START_TIMER(start) \
do { \
	(start) = MPI_Wtime();\
} while(0)

#define STOP_TIMER(phaseNo, start, end, total, rank, fptr) \
do { \
	(end) = MPI_Wtime(); \
	if ((rank) != (0)) continue; \
	fprintf((fptr), ("Phase %d completed in %f sec\n"), (phaseNo), \
	    ((end)-(start))); \
} while(0)


#define MASTER 0

int intComp(const void *, const void *);
int binarySearch(int, int, int *);
int *createIntArray(int);
int * mergePartitions(int **, int *, int);
void getLocalPartitionSizes(int *, int *, int *, int);
void partitionArray(int *, int *, int *, int *, int,  int);
void phase0(int, int, int, int *, int **);
void phase1(int, int, int, int *, int *, int **);
void phase2(int, int, int *, int **, int **);
void phase3(int, int, int, int *, int *, int **, int **, int **, int ***);
void phase4(int, int *, int **, int **);
void phase5(int, int, int, int *, int *, int **, int **);
void psrs(int, int, int, int *, char *); 
void merge(int *, int *, int *, int, int);
void tearDown();
