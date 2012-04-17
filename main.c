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

#include <time.h>
#include <pthread.h>

#include "rngs.h"
#include "rvgs.h"

#define DECLARE_TIME_MEASUREMENT struct timeval start, finish
#define START_TIME_MEASURE gettimeofday(&start, 0)
#define FINISH_TIME_MEASURE gettimeofday(&finish, 0)
#define TIME_ELAPSED ((finish.tv_sec-start.tv_sec)*1000 + (finish.tv_usec-start.tv_usec)/1000)

#define MALLOC_IT(ptr, sz) (ptr) = malloc((sz) * sizeof(*(ptr)))

/*
 * Struct that represents header of sparse matrix in
 * YSMF (CSR) format
 */
typedef struct
{
   int offset;
   int nnz;
   int m;
   int *a;
   int *ia;
   int *ja;
} csr;

// structs to be used as arguments for thread functions

typedef struct
{
	int *a_ptr;
	int *ia_ptr;
	int *ja_ptr;
	csr *mx;
	int offset;
} args_mem_copy;

typedef struct
{
	csr *ma;
	csr *mb;
	csr *out;
} args_add;


// destructors for structs

void free_csr(csr* m) {
	free(m->a);
	free(m->ia);
	free(m->ja);
}

void free_args_mem_copy(args_mem_copy* a) {
	free(a->a_ptr);
	free(a->ia_ptr);
	free(a->ja_ptr);
}

void free_args_add(args_add* a) {
	free(a->ma);
	free(a->mb);
	free(a->out);
}

#define FREE_ARRAY(ind, sz, arr, free_func) \
	int ind; \
	for (ind = 0; ind < sz; ++ind) { \
		free_func(&arr[ind]); \
	}

/*
 * Takes O(nnz) time
 */
void add(csr* ma, csr* mb, csr* out)
{
	int m = ma->m;

	int sz = ma->nnz + mb->nnz;

	MALLOC_IT(out->a, sz);
	MALLOC_IT(out->ia, m + 1);
	MALLOC_IT(out->ja, sz);

	// preprocessing is very fast ( O(1) )

	int c_iter = 0;
	int currow = 0;
	while (currow < m)
	{
		out->ia[currow] = c_iter;

		int a_from = ma->ia[currow] - ma->offset;
		int a_to = ma->ia[currow + 1] - ma->offset;
		int b_from = mb->ia[currow] - mb->offset;
		int b_to = mb->ia[currow + 1] - mb->offset;

		int a_iter = a_from;
		int b_iter = b_from;
		while ((a_iter < a_to) || (b_iter < b_to))
		{
			if ((a_iter < a_to) && (b_iter < b_to) && (ma->ja[a_iter] == mb->ja[b_iter]))
			{
				out->ja[c_iter] = ma->ja[a_iter];
				out->a[c_iter] = ma->a[a_iter] + mb->a[b_iter];
				++a_iter;
				++b_iter;
			}
			else if ((b_iter >= b_to) || ((a_iter < a_to) && (ma->ja[a_iter] < mb->ja[b_iter])))
			{
				out->ja[c_iter] = ma->ja[a_iter];
				out->a[c_iter] = ma->a[a_iter];
				++a_iter;
			}
			else
			{
				out->ja[c_iter] = mb->ja[b_iter];
				out->a[c_iter] = mb->a[b_iter];
				++b_iter;
			}
			++c_iter;
		}
		++currow;
	}

	out->m = m;
    out->ia[m] = c_iter;
    out->nnz = c_iter;
    return;
}

void *thread_func_add(void *vptr_args)
{
	args_add args = *(args_add*)vptr_args;
	add(args.ma, args.mb, args.out);
    return NULL;
}

/*
 * returns a struct with pointers to correct addresses,
 * so we get correct submatrix in place
 *
 * Takes O(1) time
 */
void get_rows(csr* m, int from_row, int to_row, csr* out)
{
	int from_index = m->ia[from_row];
	int to_index=  m->ia[to_row + 1];
	out->offset = from_index;
	out->nnz = to_index - from_index;
	out->m = to_row - from_row + 1;
	out->a = &m->a[from_index];
	out->ia = &m->ia[from_row];
	out->ja = &m->ja[from_index];
    return;
}

/**
 * Takes O(nnz+m) time
 */
void matrix_mem_copy(int* a_ptr, int* ia_ptr, int* ja_ptr, csr* mx, int offset)
{
	memcpy(a_ptr, mx->a, mx->nnz*sizeof(int));

	// memcpy with addition ( + offset )
	int  *wdst = ia_ptr;
	int  *wsrc = mx->ia;
	int i;
	for(i = 0; i < mx->m; i++)
	   *(wdst++) = *(wsrc++) + offset;

	memcpy(ja_ptr, mx->ja, mx->nnz*sizeof(int));
}

void *thread_func_mem_copy(void *vptr_args)
{
	args_mem_copy args = *(args_mem_copy*)vptr_args;
	matrix_mem_copy(args.a_ptr, args.ia_ptr, args.ja_ptr,args.mx,args.offset);
    return NULL;
}

/**
 * Takes O(1) time for preprocessing
 * and launches matrix_mem_copy,
 * which take O(nnz+m) time, but can be executed in parallel.
 */
void join_rows(size_t size, csr* matrices, csr* out)
{
	int nnz = 0;
	int m = 0;

	int i;
	for (i = 0; i < size; ++i)
	{
		csr mx = matrices[i];
		nnz += mx.nnz;
		m += mx.m;
	}

	MALLOC_IT(out->a, nnz);
	MALLOC_IT(out->ia, m + 1);
	MALLOC_IT(out->ja, nnz);

	int out_index_iter = 0;
	int out_row_iter = 0;

	pthread_t* MALLOC_IT(threads, size);
	args_mem_copy* MALLOC_IT(args, size);

	for (i = 0; i < size; ++i)
	{
		csr mx = matrices[i];

		// try done this in parallel
		args[i].a_ptr = out->a + out_index_iter;
		args[i].ia_ptr =	out->ia + out_row_iter;
		args[i].ja_ptr = out->ja + out_index_iter;
		args[i].mx = &mx;
		args[i].offset = out_index_iter;

//		cpu_set_t cpuset;
//		CPU_SET(i,&cpuset);
//		pthread_setaffinity_np(threads[i], sizeof(cpu_set_t), &cpuset);
		pthread_create(&threads[i], NULL, thread_func_mem_copy, &args[i]);

		out_index_iter += mx.nnz;
		out_row_iter += mx.m;
	}

	for (i = 0; i < size; ++i)
	{
		pthread_join(threads[i], NULL);
	}
	free(threads);
	free(args);

	out->m = m;
	out->nnz = nnz;
	out->ia[m] = nnz;
    return;
}

void mr_add(csr* ma, csr* mb, int n_threads, csr* out)
{
   DECLARE_TIME_MEASUREMENT;
   START_TIME_MEASURE;
   int m = ma->m;
   int subm = (int) ceil((m + 0.0) / (n_threads + 0.0));


   csr* MALLOC_IT(a_arr, n_threads);
   csr* MALLOC_IT(b_arr, n_threads);
   csr* MALLOC_IT(sum_arr, n_threads);

   // preprocessing takes O(1) ~ 0 ms

   pthread_t* MALLOC_IT(threads, n_threads);
   args_add* MALLOC_IT(args, n_threads);

   int from_row = 0;
   int j;
   for (j = 0; j < n_threads; ++j)
   {
	   int to_row = from_row + subm - 1;
	   if (to_row >= m)
		   to_row = m-1;

	   // map step ( O(1) ) ~ 0 ms

	   get_rows(ma, from_row, to_row, &a_arr[j]);
	   get_rows(mb, from_row, to_row, &b_arr[j]);

	   // payload ( O(nnz) ), try to do this in parallel

	   args[j].ma = &a_arr[j];
	   args[j].mb = &b_arr[j];
	   args[j].out = &sum_arr[j];

	   pthread_create(&threads[j], NULL, thread_func_add, &args[j]);
	   printf("thread %d created\n", j);
	   from_row += subm;
   }
   for (j = 0; j < n_threads; ++j)
   {
	   pthread_join(threads[j], NULL);
   }
   free(a_arr);
   free(b_arr);

   free(threads);
   free(args);
   FINISH_TIME_MEASURE;
   printf("all add()'s took %ld ms\n", TIME_ELAPSED);


   START_TIME_MEASURE;
   join_rows(n_threads, sum_arr, out);
   FINISH_TIME_MEASURE;
   printf("join_rows() took %ld ms\n", TIME_ELAPSED);

   FREE_ARRAY(ksum, n_threads, sum_arr, free_csr);
   free(sum_arr);
}

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

	out->a = (int*) malloc(nnz*sizeof(int));
	out->ia = (int*) malloc((m+1)*sizeof(int));
	out->ja = (int*) malloc(nnz*sizeof(int));

	int c_iter = 0;

	int i;
	for (i = 0; i < m; i++)
	{
		out->ia[i] = c_iter;
		int rnd_cnt;
		if (i < m - 1)

			rnd_cnt= FastBinomial(n, nnz_part);
		else
			rnd_cnt = nnz - c_iter;
		if (c_iter + rnd_cnt > nnz)
			rnd_cnt = nnz - c_iter;

		int j;
		for (j = 0; j < rnd_cnt; j++)
		{

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

void debug_print_csr(csr* m)
{
	int i;
	for (i = 0; i < m->nnz; ++i)
	{
		printf("%d ", m->a[i]);
	}
	printf("\n");

	for (i = 0; i < m->m + 1; ++i)
	{
		printf("%d ", m->ia[i]);
	}
	printf("\n");

	for (i = 0; i < m->nnz; ++i)
	{
		printf("%d ", m->ja[i]);
	}
	printf("\n\n");
}

void test() {
	int a[] = {1, 2, 3, 9, 1, 4};
	int ia[] = { 0, 2, 4, 6 };
	int ja[] = { 0, 1, 1, 2, 1, 2 };

	int a2[] = {3, 8, 1, 5, 2, 1};
	int ia2[] = { 0, 2, 4, 6 };
	int ja2[] = { 0, 2, 1, 3, 0, 1 };

	csr ma = {0, 6, 3, a, ia, ja};
	csr mb = {0, 6, 3, a2, ia2, ja2};

	debug_print_csr(&ma);
	debug_print_csr(&mb);

	csr mres;

	printf("%s\n", "GetRows");
	get_rows(&ma, 0, 1, &mres);
	debug_print_csr(&mres);


	printf("%s\n", "Addition");
	add(&ma, &mb, &mres);
	debug_print_csr(&mres);

	csr mparts[2];
	printf("%s\n", "Join");
	get_rows(&ma, 0, 1, &mparts[0]);
	get_rows(&ma, 2, 2, &mparts[1]);
	join_rows(2, mparts, &mres);
	debug_print_csr(&mres);

	printf("%s\n", "ParAdd");
	csr mparadd;
	mr_add(&ma, &mb, 2, &mparadd);
	debug_print_csr(&mparadd);

	printf("%s\n", "Generate");
	csr mgen;

	generate_csr(10000,10000,0.01,&mgen);

	return;
}


int main(int argc, char **argv)
{

	DECLARE_TIME_MEASUREMENT;

	int m = 20000, n = 20000;
	double p = 0.12;

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

//	printf("%s\n", "SeqAdd");
//	START_TIME_MEASURE;
//	add(&ma, &mb, &mout);
//	FINISH_TIME_MEASURE;
//	printf("Simple summation took %ld ms\n", TIME_ELAPSED);

	printf("%s\n", "ParAdd");
	START_TIME_MEASURE;
	mr_add(&ma, &mb, 8, &mout);
	FINISH_TIME_MEASURE;
	printf("Summation took %ld ms\n", TIME_ELAPSED);

	free_csr(&ma);
	free_csr(&mb);
	free_csr(&mout);

	return 0;
}

