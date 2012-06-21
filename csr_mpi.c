/*
 * csr_mpi.c
 *
 *  Created on: 21.06.2012
 *      Author: vbolshutkin
 */

#include <math.h>
#include "mpi.h"

#include "common.h"
#include "csr.h"
#include "csr_sum.h"
#include "csr_mpi.h"



int mpi_pack_csr_size(csr* m) {

	int size = 0, total_size = 0;

	MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &size);
	total_size += size;
	MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &size);
	total_size += size;
	MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &size);
	total_size += size;
	MPI_Pack_size(m->nnz, MPI_INT, MPI_COMM_WORLD, &size);
	total_size += size;
	MPI_Pack_size(m->m + 1, MPI_INT, MPI_COMM_WORLD, &size);
	total_size += size;
	MPI_Pack_size(m->nnz, MPI_INT, MPI_COMM_WORLD, &size);
	total_size += size;

	return total_size;
}

void mpi_pack_csr(csr* m, char* buff, int size) {

	int position = 0;

	MPI_Pack(&m->nnz, 1, MPI_INT, buff, size, &position, MPI_COMM_WORLD);
	MPI_Pack(&m->m, 1, MPI_INT, buff, size, &position, MPI_COMM_WORLD);
	MPI_Pack(&m->offset, 1, MPI_INT, buff, size, &position, MPI_COMM_WORLD);
	MPI_Pack(m->a, m->nnz, MPI_INT, buff, size, &position, MPI_COMM_WORLD);
	MPI_Pack(m->ia ,m->m + 1, MPI_INT,  buff, size, &position, MPI_COMM_WORLD);
	MPI_Pack(m->ja, m->nnz, MPI_INT,  buff, size, &position, MPI_COMM_WORLD);
}

void mpi_unpack_csr(char* buff, int size, csr* mx) {

	int position = 0;
	MPI_Unpack(buff, size, &position, &mx->nnz, 1, MPI_INT, MPI_COMM_WORLD);
	MPI_Unpack(buff, size, &position, &mx->m, 1, MPI_INT, MPI_COMM_WORLD);
	MPI_Unpack(buff, size, &position, &mx->offset, 1, MPI_INT, MPI_COMM_WORLD);

	printf("%s\n", "unpack 1/2");
	printf("%d %d %d\n", mx->nnz, mx->m, mx->offset);


	MALLOC_IT(mx->a, mx->nnz);
	MALLOC_IT(mx->ia, mx->m + 1);
	MALLOC_IT(mx->ja, mx->nnz);
	mx->subcsr_cnt = 0;

	MPI_Unpack(buff, size, &position, mx->a, mx->nnz, MPI_INT, MPI_COMM_WORLD);
	MPI_Unpack(buff, size, &position, mx->ia, mx->m + 1, MPI_INT, MPI_COMM_WORLD);
	MPI_Unpack(buff, size, &position, mx->ja, mx->nnz, MPI_INT, MPI_COMM_WORLD);
}

void mpi_csr_add_master(csr* ma, csr* mb, int n_proc, csr* out)
{

   int m = ma->m;
   int subm = (int) ceil((m + 0.0) / (n_proc + 0.0));

   csr* MALLOC_IT(a_arr, n_proc);
   csr* MALLOC_IT(b_arr, n_proc);
   csr* MALLOC_IT(sum_arr, n_proc);

   // preprocessing takes O(1) ~ 0 ms

   int from_row = 0;
   int j;
   for (j = 1; j <= n_proc; ++j) {
	   int to_row = from_row + subm - 1;
	   if (to_row >= m)
		   to_row = m-1;

	   // map step ( O(1) ) ~ 0 ms

	   get_rows(ma, from_row, to_row, &a_arr[j-1]);
	   get_rows(mb, from_row, to_row, &b_arr[j-1]);

	   // send task to other processor
	   int sz[2];
	   sz[0] = mpi_pack_csr_size(&a_arr[j-1]);
	   sz[1] = mpi_pack_csr_size(&b_arr[j-1]);

	   MPI_Send(sz, 2, MPI_INT, j, 0, MPI_COMM_WORLD);
	   printf("%s: %d, %d\n", "size sent", sz[0], sz[1]);
	   char* MALLOC_IT(outbuff, sz[0]);
	   mpi_pack_csr(&a_arr[j-1], outbuff, sz[0]);
	   printf("%s\n", "pack done");
	   MPI_Send(outbuff, sz[0], MPI_PACKED, j, 1, MPI_COMM_WORLD);
	   printf("%s\n", "send done");

	   REALLOC_IT(outbuff, sz[1]);
	   mpi_pack_csr(&b_arr[j-1], outbuff, sz[1]);
	   printf("%s\n", "pack2 done");
	   MPI_Send(outbuff, sz[1], MPI_PACKED, j, 2, MPI_COMM_WORLD);
	   printf("%s\n", "send2 done");

	   free(outbuff);

	   from_row += subm;
   }
   for (j = 1; j <= n_proc; ++j) {
	   printf("%s\n", "start sync");
	   MPI_Status stat;
	   int insz = 0;
	   char* inbuff;
	   MPI_Recv(&insz, 1, MPI_INT, j, 3, MPI_COMM_WORLD, &stat);
	   printf("master sync: size %d\n", insz);
	   MALLOC_IT(inbuff, insz);
	   MPI_Recv(inbuff, insz, MPI_PACKED, j, 4, MPI_COMM_WORLD, &stat);
	   printf("master sync: got answer\n");
	   mpi_unpack_csr(inbuff, insz, &sum_arr[j-1]);
	   printf("master sync: unpack\n");

	   free(inbuff);
	   printf("master sync: free\n");
   }
   free(a_arr);
   free(b_arr);

   printf("master after sync: before join\n");
   join_rows(n_proc, sum_arr, out);
}

void mpi_csr_add_slave()
{
	int size[2];
	char* buff, *outbuff;

	MPI_Status stat;

	printf("slave: %s\n", "start work");

	MPI_Recv(size, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
	printf("%s: %d, %d\n", "size received", size[0], size[1]);

	MALLOC_IT(buff, size[0]);

	csr a, b, out;

	MPI_Recv(buff, size[0], MPI_PACKED, 0, 1, MPI_COMM_WORLD, &stat);
	mpi_unpack_csr(buff, size[0], &a);

	REALLOC_IT(buff, size[1]);

	MPI_Recv(buff, size[1], MPI_PACKED, 0, 2, MPI_COMM_WORLD, &stat);
	mpi_unpack_csr(buff, size[1], &b);

	free(buff);

	printf("slave: %s\n", "before add");
	add(&a, &b, &out);
	printf("slave: %s\n", "after add");

	printf("slave: %s %d\n", "a offset = ", a.offset);
	printf("slave: %s %d\n", "b offset = ", b.offset);
	printf("slave: %s %d\n", "out offset = ", out.offset);

	int sz = mpi_pack_csr_size(&out);
	printf("slave: %s %d\n", "response size", sz);
	MPI_Send(&sz, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);

	MALLOC_IT(outbuff, sz);
	mpi_pack_csr(&out, outbuff, sz);

	MPI_Send(outbuff, sz, MPI_PACKED, 0, 4, MPI_COMM_WORLD);
	printf("slave: %s\n", "response sent");
	free(outbuff);
}
