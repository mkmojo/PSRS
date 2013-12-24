/*
 * Implementations of the PSRS algorithm using the Message Passing Interface.
 * Marko Babic - mbabic - Cmput 481/681
 */
#include "psrs.h"

static int validateResults;

void
psrs(int size, int rank, int nInts, int *toSort, char *fname)
{
	char hname[256];
	double start, end, total = 0, algStart;
	int i, *privateInts, *regularSamples, *collectedSamples, *pivots, 
	    *partitionIndices, *localPartitionSizes, *incomingPartitionSizes, 
	    **partitions, *mergedPartitions, *partitionSizes, *sortedArray; 

	FILE *fptr = NULL;
	if (rank == MASTER) fptr = fopen(fname, "a");
	
	memset (hname, '\0', sizeof(unsigned char)*256);
	gethostname(hname, 255);

	DBPRINT(("%d of %d running on pid %d %s\n", rank, size, getpid(), 
	    hname));

	if (rank == MASTER) {
		collectedSamples = calloc(1, sizeof(int)*size*size);
		if (!collectedSamples) {
			err(1,"Failed to allocate memory for collected samples "
			   "array for master process");
			MPI_Finalize();
			exit(1);
		}				
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if (nInts < VALIDATION_THRESHOLD) validateResults = 1;

	/* 
	 * "Phase 0" in which the array is split into contiguous chunks
	 * distributed amongst the processes.
	 */
	phase0(size, rank, nInts, toSort, &privateInts);

	/* Phase 1 */

	algStart = MPI_Wtime();
	START_TIMER((start));
	phase1(size, rank, nInts, toSort, privateInts, &regularSamples);
	MPI_Barrier(MPI_COMM_WORLD);
	STOP_TIMER((1), (start), (end), (total), (rank), (fptr));

	/* Phase 2 */
	START_TIMER((start));
	phase2(size, rank, regularSamples, &collectedSamples, &pivots);
	MPI_Barrier(MPI_COMM_WORLD);
	STOP_TIMER((2), (start), (end), (total), (rank), (fptr));

	/* Phase 3 */
	START_TIMER((start));
	phase3(size, rank, nInts, privateInts, pivots, &localPartitionSizes,
	    &partitionIndices, &incomingPartitionSizes, &partitions); 
	MPI_Barrier(MPI_COMM_WORLD);
	STOP_TIMER((3), (start), (end), (total), (rank), (fptr));

	/* Phase 4 */
	START_TIMER((start));
	phase4(size, incomingPartitionSizes, partitions, &mergedPartitions);
	MPI_Barrier(MPI_COMM_WORLD);
	STOP_TIMER((4), (start), (end), (total), (rank), (fptr));
	if (rank == MASTER) {
		total = end - algStart;
		fprintf(fptr, "The algorithm took %f seconds in total\n", total);
	}

	if (!validateResults) return;


	/* Run validation test */
	
	/* "Phase 5" concatenate the lists back at the master */
	phase5(size, rank, nInts, incomingPartitionSizes,
	    mergedPartitions, &partitionSizes, &sortedArray);

	/* 
	 * Assert that the array is equivalent to the sorted original array
	 * where we sort the original array using a known, proven, sequential,
	 * method.
	 */
	if (rank == MASTER) {
		qsort(toSort, nInts, sizeof(int), intComp);
	}

	for (i = 0; i < nInts; i++) {
		if (rank == MASTER) {
			if (toSort[i] != sortedArray[i]) {
				printf("OH NO, got %d at pos %d, expected "
				    "%d\n", sortedArray[i], i, toSort[i]);
			}
		}
	}
	
	if (rank == MASTER) fclose(fptr);
}


/*
 * Calculates the sizes of the partitions of the local subarray based on
 * the indices at which each partition begins/ends.
 */
void
getLocalPartitionSizes(int *localPartitionSizes, int *partitionIndices,
    int * privateInts, int size)
{
	int partitionSize, i;

	for (i = 0; i < size; i++) {
		if (partitionIndices[i] == -1) partitionSize = 0;
		else if (i == 0) partitionSize = partitionIndices[i] + 1;
		else {
			partitionSize = partitionIndices[i] -
			    partitionIndices[i-1];
		}
		localPartitionSizes[i] = partitionSize;
	}
}

/*
 * Populates the subarraySizes array based on the position of the pivots
 * in the given privateInts array.
 * partitionIndices works as follows:
 * partitionIndices[i] = d => every element to the left of privateInts[d]
 * (including privateInts[d] itself) down to but NOT including element 
 * privateInts[partitionIndices[i - 1]].  If d == -1 or (d - partitionIndices[
 * i-1]) == 0, then the partition is empty. 
 */
void
partitionArray(int *privateInts, int *pivots, int *partitionIndices, 
    int * localPartitionSizes, int nInts, int size)
{
	int i, key, index;
	partitionIndices[size - 1] = (nInts / size) - 1;

	for (i = 0; i < size - 1; i++) {
		key = pivots[i];
		index = binarySearch(key, nInts / size, privateInts);
		partitionIndices[i] = index;
	}

	getLocalPartitionSizes(localPartitionSizes, partitionIndices,
	    privateInts, size);	
}


/*
 * Integer comparason functions used for call to qsort() in phase 1 of the 
 * PSRS algorithm.
 */
int
intComp(const void *a, const void *b)
{
	int *x = (int *)a, *y = (int *)b;
	if (*x > *y) return 1;
	if (*x < *y) return -1;
	return 0;
}

/*
 * Binary search for index in array such that array[index] == key or 
 * array[index] is element of greatest index with value < key.
 */
int
binarySearch(int key, int size, int *array)
{
	int low = 0, high = size - 1, mid;
	if (array[0] > key) return -1;

	mid = (low + high) / 2;
	while (low <= high) {
		if (array[mid] > key) high = mid - 1;
		else if (array[mid] < key) low = mid + 1;
		else break;
		mid = (low + high) / 2;
	}	

	if (mid == size - 1) return mid;
	while (array[mid] == array[mid + 1]) {
		mid++;
		if (mid == size-1) return mid;
	}
	return mid;
}

/*
 * Performs a sorting merge on the given array of partitions and returns
 * the result.
 */
int *
mergePartitions(int **partitions, int *incomingPartitionSizes, int size)
{
	int *ret = NULL, *retCpy = NULL, mergedSize = 0, oldMergedSize = 0, i;

	ret = partitions[0];
	mergedSize += incomingPartitionSizes[0];
	for (i = 1; i < size; i++) {
		oldMergedSize = mergedSize;
		mergedSize += incomingPartitionSizes[i];
		retCpy = ret;
		ret = calloc(mergedSize, sizeof(int));
		merge(partitions[i], retCpy, ret, incomingPartitionSizes[i],
		    oldMergedSize);
	}
	return ret;
}

void
merge(int *a, int *b, int *ret, int aLen, int bLen) 
{
	int i = 0, k = 0, j = 0;

	for (i = 0; i < aLen + bLen; i++) {
		if ( (j < aLen) && (k < bLen) ) {
			if (a[j] < b[k]) {
				ret[i] = a[j];
				j++;
			} else {
				ret[i] = b[k];
				k++;
			}
		} else if (j == aLen) {
			ret[i] = b[k];
			k++;
		} else {
			ret[i] = a[j];
			j++;
		}
	}

}


void
phase0(int size, int rank, int nInts, int *toSort, int **privateInts)
{
	*privateInts = calloc((nInts / size), sizeof(int));
	if (!*privateInts) {
		err(1, "Proc %d could not allocate memory for private subarray"
		    , rank);
		MPI_Finalize();
		exit(1);
	}

	/* Master scatters partitions to each process for local sorting */
	MPI_Scatter(toSort, nInts / size, MPI_INT, *privateInts, nInts / size, 
	    MPI_INT, MASTER, MPI_COMM_WORLD);
}


/*
 * A subpartition of the array to be sorted is scattered to each node, the 
 * node then sorts its partition and selects regular samples.
 */
void
phase1(int size, int rank, int nInts, int *toSort, int *privateInts,
    int **regularSamples)
{
	int w, i;
	/* Allocate memory for local phase1 arrays */
	*regularSamples = calloc(size, sizeof(int));		/* Process sorts local subarray */
	qsort(privateInts, nInts / size, sizeof(int), intComp);

#ifdef DEBUG1
	/* make sure local subarray is actually sorted */
	for (i = 0; i < nInts / size - 1; i++) {
		if (privateInts[i] > privateInts[i+1]) {
			printf("\n %d's subarray not sorted!", rank);
			break;
		}
	}
#endif

	w = nInts / (size * size);	
	/* Each process selects regular samples of its subarray */
	for (i = 0; i < size; i++) {
		(*regularSamples)[i] = privateInts[i*w];
	}
} /* Phase 1 Complete */



void
phase2(int size, int rank, int *regularSamples, int **collectedSamples,
    int **pivots)
{
	int rho, i;

	if (rank == MASTER) *collectedSamples = calloc(size*size, sizeof(int));
	*pivots = calloc(size, sizeof(int));

	MPI_Gather(regularSamples, size, MPI_INT, *collectedSamples, size,
	    MPI_INT, MASTER, MPI_COMM_WORLD);

	free(regularSamples);
#ifdef DEBUG1
	if (rank == MASTER) {
		for (i = 0; i < size*size; i++) {
				printf("%d : %d\n", i, (*collectedSamples)[i]);
			}
	}
#endif

	/* Master processor sorts the regular samples */
	if (rank == MASTER) {
		qsort(*collectedSamples, size * size, sizeof(int), intComp);
	}

	rho = size / 2;

	if (rank == MASTER) {
		for(i = 0; i < size - 1; i++) {
			(*pivots)[i] = (*collectedSamples)[(size*(i+1)) + rho];
		}
	}

	/* Copy of pivots sent to each process */
	MPI_Bcast(*pivots, size - 1, MPI_INT, MASTER, MPI_COMM_WORLD);	

} /* Phase 2 complete */



void
phase3(int size, int rank, int nInts, int *privateInts, int *pivots,
    int **localPartitionSizes, int **partitionIndices, 
    int **incomingPartitionSizes, int ***partitions)
{
	MPI_Request *sendRequests, *recRequests;
	MPI_Status *sendStatuses, *recStatuses;

	int i;

	sendRequests = calloc(size, sizeof(MPI_Request));
	recRequests = calloc(size, sizeof(MPI_Request));
	sendStatuses = calloc(size, sizeof(MPI_Status));
	recStatuses = calloc(size, sizeof(MPI_Status));
	
	*localPartitionSizes = calloc(size, sizeof(int));
	*partitionIndices = calloc(size, sizeof(int));
	*incomingPartitionSizes = calloc(size, sizeof(int));
	*partitions = calloc(size, sizeof(int *));
	
	partitionArray(privateInts, pivots, *partitionIndices,
	    *localPartitionSizes, nInts, size);
#ifdef HELLO 
	printf("Hello, I am %d.  This is my privateInts array:\n", rank);
	for (i = 0; i < nInts / size; i++) {
		printf("%d: privateInts[%d] = %d\n", rank, i, privateInts[i]);
	}
	printf("Hello, I am %d.  These are the pivots I have:\n", rank);
	for (i = 0; i < size - 1; i++) {
		printf("%d: pivots[%d] = %d\n", rank, i, pivots[i]);
	}
	printf("Hello, I am %d.  These are my partition indices:\n", rank);
	for (i = 0; i < size; i++) {
		printf("%d: partitionIndices[%d] = %d\n", rank, i,
		    partitionIndices[i]);
	}	
	printf("Hello, I am %d.  These are my partition sizes:\n", rank);
	for (i = 0; i < size; i++) {
		printf("%d: localPartitionSizes[%d] = %d\n", rank, i,
		    localPartitionSizes[i]);
	}	
#endif

	/* 
	 * Each process shares with the other processes the size of the 
	 * partition it will be sending that process so that the receiving
	 * process can allocate buffers of appropriate sizes.
	 */
	MPI_Alltoall(*localPartitionSizes, 1, MPI_INT, *incomingPartitionSizes,
	    1, MPI_INT, MPI_COMM_WORLD);
	/* 
	 * Now that each process knows what size partition to expect, it can
	 * appropriately allocate buffers.
	 */
	for (i = 0; i < size; i++) {
		(*partitions)[i] = calloc((*incomingPartitionSizes)[i], 
		    sizeof(int));
	}
	
	/* We are now ready for the processes to exchange partitions. */
	for (i = 0; i < size; i++) {
		
		MPI_Irecv((*partitions)[i],(*incomingPartitionSizes)[i],MPI_INT,
		    i, 0, MPI_COMM_WORLD, &(recRequests[i]));
		
		if (i == 0) {
			MPI_Isend(privateInts,(*localPartitionSizes)[i],MPI_INT, 			    i, 0, MPI_COMM_WORLD, &(sendRequests[i]));
		} else {
			MPI_Isend(privateInts + (*partitionIndices)[i-1] + 1,
			    (*localPartitionSizes)[i], MPI_INT, i, 0,
			    MPI_COMM_WORLD, &(sendRequests[i]));
		}
	}
	MPI_Waitall(size, recRequests, recStatuses);
	MPI_Waitall(size, sendRequests, sendStatuses);


	free(*localPartitionSizes);
	free(privateInts);
	free(*partitionIndices);
	free(recRequests);
	free(sendRequests);
	free(recStatuses);
	free(sendStatuses);

#ifdef DEBUG1
	int j;
	for (i = 0; i < size; i++) {
		for (j = 0; j < incomingPartitionSizes[i]; j++) {
			printf("%d: partition[%d][%d] = %d\n", rank, i, j,
			    partitions[i][j]);
		}
	}
#endif
} /* Phase 3 complete */

void
phase4(int size, int *incomingPartitionSizes, 
    int **partitions, int **mergedPartitions)
{
	*mergedPartitions = mergePartitions(partitions, 
	    incomingPartitionSizes, size);
	free(*partitions);
}

void
phase5(int size, int rank, int nInts, int *incomingPartitionSizes,
    int *mergedPartitions, int **partitionSizes, int ** sortedArray)
{
	MPI_Request sendRequest, *recRequests;
	MPI_Status sendStatus, *recStatuses;
	int mergedSize = 0, nextPos = 0, i;

	recRequests = calloc(size, sizeof(MPI_Request));
	recStatuses = calloc(size, sizeof(MPI_Status));
	for (i = 0; i < size; i++) {
		mergedSize += incomingPartitionSizes[i];
	}


	if (rank == MASTER) {
		*sortedArray = calloc(nInts, sizeof(int));
		*partitionSizes = calloc(size, sizeof(int));
	}
	
	MPI_Gather(&mergedSize, 1, MPI_INT, *partitionSizes, 1, MPI_INT, MASTER,
	    MPI_COMM_WORLD);

	if (rank == MASTER) {
		for (i = 0; i < size; i++) {
			if (i == 0) MPI_Irecv(*sortedArray, 
			    (*partitionSizes)[i], MPI_INT, i, 0, 
			    MPI_COMM_WORLD, &(recRequests[i]));
			else MPI_Irecv((*sortedArray) + nextPos,
			    (*partitionSizes)[i],MPI_INT, i, 0, MPI_COMM_WORLD, 
			    &(recRequests[i]));
			nextPos += (*partitionSizes)[i];
		}
	}	
	MPI_Isend(mergedPartitions, mergedSize, MPI_INT, MASTER, 0, 
	    MPI_COMM_WORLD, &sendRequest);
	if (rank == MASTER) MPI_Waitall(size, recRequests, recStatuses);
	MPI_Wait(&sendRequest, &sendStatus);

	free(mergedPartitions);
	free(incomingPartitionSizes);
	free(recRequests);
	free(recStatuses);
}
