/*
 * main.c
 *
 *  Created on: 15.04.2012
 *      Author: vbolshutkin
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include <sys/time.h>
#include <pthread.h>

#include <unistd.h>

#include "mpi.h"

#include "common.h"


#include "rngs.h"
#include "rvgs.h"

#include "csr.h"
#include "csr_sum.h"
#include "csr_mpi.h"

#define DECLARE_TIME_MEASUREMENT struct timeval start, finish
#define START_TIME_MEASURE gettimeofday(&start, 0)
#define FINISH_TIME_MEASURE gettimeofday(&finish, 0)
#define TIME_ELAPSED ((finish.tv_sec-start.tv_sec)*1000 + (finish.tv_usec-start.tv_usec)/1000)

/*
 * fast (O(1)) but approximate implementation
 * Binomial worked in O(n) -- too long
 */
long FastBinomial(long n, double p)
{
  return round(Normal(n * p, n * p * (1 - p)));
}


void generate_csr(long n, long m, float nnz_part, csr* out)
{
	PutSeed(time(NULL));

	long nnz = ceil(n*m*nnz_part);

	out->subcsr_cnt = 0;
	out->subcsrs = NULL;

	out->a = (int*) malloc(nnz*sizeof(int));
	out->ia = (int*) malloc((m+1)*sizeof(int));
	out->ja = (int*) malloc(nnz*sizeof(int));

	int c_iter = 0;

	int i;
	for (i = 0; i < m; i++) {
		out->ia[i] = c_iter;
		int rnd_cnt;
		if (i < m - 1)

			rnd_cnt= FastBinomial(n, nnz_part);
		else
			rnd_cnt = nnz - c_iter;
		if (c_iter + rnd_cnt > nnz)
			rnd_cnt = nnz - c_iter;

		int j;
		for (j = 0; j < rnd_cnt; j++) {

			out->ja[c_iter] = Equilikely(j*(n-1)/rnd_cnt, (j+1)*(n-1)/rnd_cnt);
			out->a[c_iter] = Equilikely(1, INT_MAX / 2);
			++c_iter;
		}
	}
	out->ia[m] = c_iter;
	out->nnz = c_iter;
	out->m = m;
	out->offset = 0;
	return;
}




int main(int argc, char **argv)
{
	int m = 20000, n = 20000;
	double p = 0.10;

	int c;
	while ((c = getopt (argc, argv, "hvn:m:p:")) != -1) {
		 switch (c)
		   {
		   case 'h':
			 break;
		   case 'v':
			 break;
		   case 'n':
			 sscanf(optarg, "%d", &n);
			 if (n < 1 || n > INT_MAX) {
				 return -1;
			 }
			 break;
		   case 'm':
			 sscanf(optarg, "%d", &m);
			 if (m < 1 || m > INT_MAX) {
				 return -1;
			 }
			 break;
		   case 'p':
			 sscanf(optarg, "%lf", &p);
			 if (p <= 0 || p >= 1) {
				 return -1;
			 }
			 break;
		   default:
			 abort ();
			 break;
		   }
	}


	int myid;
	int n_proc;

	/* MPI programs start with MPI_Init; all 'N' processes exist thereafter */
	MPI_Init(&argc,&argv);
	/* find out how big the SPMD world is */
	MPI_Comm_size(MPI_COMM_WORLD,&n_proc);
	/* and this processes' rank is */
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	if (myid == 0) {

		DECLARE_TIME_MEASUREMENT;

		csr ma, mb, mout;

		printf("%s\n", "Generation");
		START_TIME_MEASURE;
		generate_csr(n,m,p,&ma);
		FINISH_TIME_MEASURE;
		printf("Matrix A generation took %ld ms\n", TIME_ELAPSED);

		START_TIME_MEASURE;
		generate_csr(n,m,p,&mb);
		FINISH_TIME_MEASURE;
		printf("Matrix B generation took %ld ms\n", TIME_ELAPSED);

		printf("%s\n", "SeqAdd");
		START_TIME_MEASURE;
		add(&ma, &mb, &mout);
		FINISH_TIME_MEASURE;
		double simple_time = TIME_ELAPSED;
		printf("Simple summation took %ld ms\n", TIME_ELAPSED);

		printf("%s\n", "Scatter");
		START_TIME_MEASURE;
		mpi_csr_scatter(&ma, &mb, n_proc - 1); // number of slaves
		FINISH_TIME_MEASURE;
		printf("Scattering took %ld ms\n", TIME_ELAPSED);

		printf("%s\n", "ParAdd");
		START_TIME_MEASURE;
		mpi_csr_add(n_proc - 1);
		FINISH_TIME_MEASURE;
		double mpi_add_time = TIME_ELAPSED;
		printf("Summation took %ld ms\n", TIME_ELAPSED);

		printf("With %d threads it is %.2f times faster\n", n_proc, (simple_time + 0.0) / (mpi_add_time + 0.0));

		mpi_csr_free(n_proc - 1);

		free_csr(&ma);
		free_csr(&mb);
		free_csr(&mout);
	} else {
		mpi_csr_add_slave();
	}

	MPI_Finalize();
	return 0;
}

