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

// destructor for args_mem_copy struct
void free_args_mem_copy(args_mem_copy* a)
{
	free(a->a_ptr);
	free(a->ia_ptr);
	free(a->ja_ptr);
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

/*
 * delegates call to matrix_mem_copy()
 */
void *thread_func_mem_copy(void *vptr_args)
{
	args_mem_copy args = *(args_mem_copy*)vptr_args;
	matrix_mem_copy(args.a_ptr, args.ia_ptr, args.ja_ptr,args.mx,args.offset);
    return NULL;
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

	pthread_t* MALLOC_IT(threads, in->subcsr_cnt);
	args_mem_copy* MALLOC_IT(args, in->subcsr_cnt);

	for (i = 0; i < in->subcsr_cnt; ++i)
	{
			csr mx = in->subcsrs[i];

			// try done this in parallel
			args[i].a_ptr = out->a + out_index_iter;
			args[i].ia_ptr =        out->ia + out_row_iter;
			args[i].ja_ptr = out->ja + out_index_iter;
			args[i].mx = &mx;
			args[i].offset = out_index_iter;

			pthread_create(&threads[i], NULL, thread_func_mem_copy, &args[i]);

			out_index_iter += mx.nnz;
			out_row_iter += mx.m;
	}

	for (i = 0; i < in->subcsr_cnt; ++i)
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
