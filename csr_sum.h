#ifndef CSR_SUM_H
#define CSR_SUM_H

#include "csr.h"

typedef struct {
	csr *ma;
	csr *mb;
	csr *out;
} args_add;

void free_args_add(args_add* a);

void add(csr* ma, csr* mb, csr* out);

void *thread_func_add(void *vptr_args);

void get_rows(csr* m, int from_row, int to_row, csr* out);

void join_rows(size_t size, csr* matrices, csr* out);

void mr_add(csr* ma, csr* mb, int n_threads, csr* out);

#endif
