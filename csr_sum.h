#ifndef CSR_SUM_H
#define CSR_SUM_H

#include "csr.h"

void add(csr* ma, csr* mb, csr* out);

void get_rows(csr* m, int from_row, int to_row, csr* out);

void join_rows(size_t size, csr* matrices, csr* out);

void mr_add(csr* ma, csr* mb, int n_threads, csr* out);

#endif
