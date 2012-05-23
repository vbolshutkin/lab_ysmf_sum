#ifndef CSR_H
#define CSR_H

#include "common.h"

/*
 * Struct that represents header of sparse matrix in
 * YSMF (CSR) format
 */
typedef struct csrT {
   int nnz;
   int m;

   // only for packed csr
   int offset;
   int *a;
   int *ia;
   int *ja;

   // only for fragmented csr
   int subcsr_cnt;
   struct csrT *subcsrs;
} csr;

// getters for fragmented csr. work in O(nFragments)
int get_csr_a(csr* m, int index);
int get_csr_ia(csr* m, int index);
int get_csr_ja(csr* m, int index);

// destructor
void free_csr(csr* m);

/**
 * Takes O(nnz+m) time
 */
void matrix_mem_copy(int* a_ptr, int* ia_ptr, int* ja_ptr, csr* mx, int offset);

/**
 * Function, that packs a fragmented csr
 *
 * Takes O(nFragments) time for preprocessing
 * and launches matrix_mem_copy,
 * which take O(nnz+m) time, but can be executed in parallel.
 */
void pack_csr(csr* in, int n_threads, csr* out);

#endif
