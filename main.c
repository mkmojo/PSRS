/*
 * Test runner for PSRS timing tests.
 * CMput 481/681 - Marko Tomislav Babic - mbabic
 */

#include "psrs.h"
#include "sequential_sort.h"

#define NUM_RUNS 3 
#define RSEED 92349278	/* seed for the prng */

static int *createTestArray(int); 
static void runTests(int, char **);

int
main(int argc, char **argv)
{
	runTests(argc, argv);
	return 0;
}

static void
runTests(int argc, char **argv)
{
	int *toSort = NULL;
	FILE *fptr = NULL;
	int i, size, rank, nInts;
	

	for (i = 0; i < NUM_RUNS; i++) {
		if (i == 0) {
			MPI_Init(&argc, &argv);
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			MPI_Comm_size(MPI_COMM_WORLD, &size);
		}

		if (argc != 3) {
			if (rank == MASTER) fprintf(stderr, "Usage: ./psrsTest "
			    " nInts outputfile\n");
			MPI_Finalize();
			exit(1);
		}

		nInts = atoi(argv[1]);

		if (rank == MASTER) fptr = fopen(argv[2], "a");


		if (!nInts) {
			err(1, "BOOM, your 0 integer array just got sorted.");
		}

		if (nInts % size != 0) {
			if (rank == MASTER) err(1, "Please select a number of"
			    " integers which is divisible by the number of "
			    "processes (for analytical simplicity)");
			MPI_Finalize();
			exit(1); 
		}	

		if (nInts % (size*size) != 0) {
			if (rank == MASTER) err(1, "Please select a number of "
			    "integers which is divisible by the number of "
			    "processes^2 (for analytical simplicity)");
			MPI_Finalize();
			exit(1); 
		}

		if (rank == MASTER) toSort = createTestArray(nInts);
		MPI_Barrier(MPI_COMM_WORLD);
		if (rank == MASTER) {
			fprintf(fptr, 
			    "- Test Start -------------------------------\n");
			fprintf(fptr, "Sorting %d integers using %d "
			    "processes\n", nInts, size);
			fclose(fptr);
		}

		if (size == 1) sequential_sort(nInts, toSort, argv[2]);
		else psrs(size, rank, nInts, toSort, argv[2]);	
		
		if (rank == MASTER) {
			fptr = fopen(argv[2], "a");	
		    	fprintf(fptr, 
		    	    "- Test End ---------------------------------\n");
			fclose(fptr);
		}

	}
	MPI_Finalize();	
}

/*
 * Creates and return a pointer to an array of nInts random integers.
 */
static int *
createTestArray(int nInts) 
{
	int *ints = NULL, i;
	srandom(RSEED);

	ints = calloc(1, sizeof(int)*nInts); 
	if (!ints) {
		err(1, "Failed to allocate memory for list of integers to"
		    " sort!");
		exit(EXIT_FAILURE);
	}


	for (i = 0; i < nInts; i++) {
		ints[i] = random();
	}

	return ints;
}
