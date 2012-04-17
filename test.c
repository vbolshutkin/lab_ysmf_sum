#include <stdio.h>

#include "test.h"

void debug_print_csr(csr* m)
{
	int i;
	for (i = 0; i < m->nnz; ++i) {
		printf("%d ", get_csr_a(m, i));
	}
	printf("\n");

	for (i = 0; i < m->m + 1; ++i) {
		printf("%d ", get_csr_ia(m, i));
	}
	printf("\n");

	for (i = 0; i < m->nnz; ++i) {
		printf("%d ", get_csr_ja(m, i));
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

	csr ma = {6, 3, 0, a, ia, ja, 0, NULL};
	csr mb = {6, 3, 0, a2, ia2, ja2, 0, NULL};

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
	generate_csr(10,10,0.1,&mgen);
	debug_print_csr(&mgen);

	printf("%s\n", "Test OK");

	return;
}
