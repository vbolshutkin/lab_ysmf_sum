#include <stdlib.h>
#include <string.h>
#include "csr.h"

// getters for fragmented csr. work in O(nFragments)

int get_csr_a(csr* m, int index)
{
	if (m->subcsr_cnt == 0)
		return m->a[index];

	csr* subcsr = m->subcsrs;
	while (index > subcsr->nnz) {
		index -= subcsr->nnz;
		++subcsr;
	}
	return subcsr->a[index];
}

int get_csr_ia(csr* m, int index)
{
	if (m->subcsr_cnt == 0)
		return m->ia[index];

	int global_offset = 0;
	csr* subcsr = m->subcsrs;
	while (index > subcsr->m) {
		index -= subcsr->m;
		global_offset += subcsr->nnz;
		++subcsr;
	}
	return subcsr->ia[index] + global_offset;
}

int get_csr_ja(csr* m, int index)
{
	if (m->subcsr_cnt == 0)
		return m->ja[index];

	csr* subcsr = m->subcsrs;
	while (index > subcsr->nnz) {
		index -= subcsr->nnz;
		++subcsr;
	}
	return subcsr->ja[index];
}

// csr destructor
void free_csr(csr* m)
{
	if (m->subcsr_cnt == 0) {
		free(m->a);
		free(m->ia);
		free(m->ja);
	}
	int i;
	for (i = 0; i < m->subcsr_cnt; ++i) {
		free_csr(& m->subcsrs[i]);
	}
	free(m->subcsrs);
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

/**
 * Function, that packs a fragmented csr
 *
 * Takes O(nFragments) time for preprocessing
 * and launches matrix_mem_copy,
 * which take O(nnz+m) time, but can be executed in parallel.
 */
void pack_csr(csr* in, int n_threads, csr* out)
{
	int nnz = 0;
	int m = 0;

	int i;
	for (i = 0; i < in->subcsr_cnt; ++i)
	{
			csr mx = in->subcsrs[i];
			nnz += mx.nnz;
			m += mx.m;
	}

	MALLOC_IT(out->a, nnz);
	MALLOC_IT(out->ia, m + 1);
	MALLOC_IT(out->ja, nnz);

	int out_index_iter = 0;
	int out_row_iter = 0;

#pragma omp parallel for
	for (i = 0; i < in->subcsr_cnt; ++i)
	{
			csr mx = in->subcsrs[i];

			matrix_mem_copy(out->a + out_index_iter,
					out->ia + out_row_iter,
					out->ja + out_index_iter,
					&mx,
					out_index_iter);

			out_index_iter += mx.nnz;
			out_row_iter += mx.m;
	}

	out->m = m;
	out->nnz = nnz;
	out->ia[m] = nnz;

    return;
}
