#include <err.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "mpi.h"

static int intComp(const void *, const void *);

void
sequential_sort(int nInts, int *toSort, char *fname)
{
	double start, end;
	FILE *fp = fopen(fname, "a");
	if (!fp) err(1, "Failed to open file %s", fname);

	start = MPI_Wtime();
	qsort(toSort, nInts, sizeof(int), intComp);
	end = MPI_Wtime();
	fprintf(fp, "Sorting %d ints sequentially took %f seconds.\n", nInts,
	    end-start);
	fclose(fp);
}

static int
intComp(const void *a, const void *b)
{
	int *x = (int *)a, *y = (int *)b;
	if (*x > *y) return 1;
	if (*x < *y) return -1;
	return 0;
}
