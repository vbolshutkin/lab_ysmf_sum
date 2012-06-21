#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <omp.h>

#include "csr_sum.h"

/*
 * Takes O(nnz) time. Works only with packed structs
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
	while (currow < m) {
		out->ia[currow] = c_iter;

		int a_from = ma->ia[currow] - ma->offset;
		int a_to = ma->ia[currow + 1] - ma->offset;
		int b_from = mb->ia[currow] - mb->offset;
		int b_to = mb->ia[currow + 1] - mb->offset;

		int a_iter = a_from;
		int b_iter = b_from;
		while ((a_iter < a_to) || (b_iter < b_to)) {
			if ((a_iter < a_to) && (b_iter < b_to) && (ma->ja[a_iter] == mb->ja[b_iter])) {
				out->ja[c_iter] = ma->ja[a_iter];
				out->a[c_iter] = ma->a[a_iter] + mb->a[b_iter];
				++a_iter;
				++b_iter;
			} else if ((b_iter >= b_to) || ((a_iter < a_to) && (ma->ja[a_iter] < mb->ja[b_iter]))) {
				out->ja[c_iter] = ma->ja[a_iter];
				out->a[c_iter] = ma->a[a_iter];
				++a_iter;
			} else {
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
    out->offset = 0;

	REALLOC_IT(out->a, c_iter);
	REALLOC_IT(out->ja, c_iter);

    out->subcsr_cnt = 0;
    out->subcsrs = NULL;
    return;
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
	out->subcsr_cnt = 0;
	out->subcsrs = NULL;
    return;
}

// takes O(nCPU)
void join_rows(size_t size, csr* matrices, csr* out)
{
	out->subcsr_cnt = size;
	out->subcsrs = matrices;

	int nnz = 0;
	int m = 0;

	int i;
	for (i = 0; i < size; ++i) {
		csr mx = matrices[i];
		nnz += mx.nnz;
		m += mx.m;
	}
	out->m = m;
	out->nnz = nnz;
	out->offset = 0;
}

void mr_add(csr* ma, csr* mb, int n_threads, csr* out)
{
	omp_set_dynamic(0);
	omp_set_num_threads(n_threads);

   int m = ma->m;
   int subm = (int) ceil((m + 0.0) / (n_threads + 0.0));

   csr* MALLOC_IT(a_arr, n_threads);
   csr* MALLOC_IT(b_arr, n_threads);
   csr* MALLOC_IT(sum_arr, n_threads);

   // preprocessing takes O(1) ~ 0 ms

   int from_row = 0;
   int j;
#pragma omp parallel for
   for (j = 0; j < n_threads; ++j) {
	   int to_row = from_row + subm - 1;
	   if (to_row >= m)
		   to_row = m-1;

	   // map step ( O(1) ) ~ 0 ms

	   get_rows(ma, from_row, to_row, &a_arr[j]);
	   get_rows(mb, from_row, to_row, &b_arr[j]);

	   // payload ( O(nnz) ), try to do this in parallel

	   add(&a_arr[j], &b_arr[j], &sum_arr[j]);

	   from_row += subm;
   }
   free(a_arr);
   free(b_arr);

   join_rows(n_threads, sum_arr, out);
}
