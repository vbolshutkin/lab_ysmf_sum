/*
 * csr_mpi.h
 *
 *  Created on: 21.06.2012
 *      Author: vbolshutkin
 */

#ifndef CSR_MPI_H_
#define CSR_MPI_H_

int mpi_pack_csr_size(csr* m);
void mpi_pack_csr(csr* m, char* buff, int size);
void mpi_unpack_csr(char* buff, int size, csr* mx);

void mpi_csr_add_master(csr* ma, csr* mb, int n_proc, csr* out);
void mpi_csr_add_slave();

#endif /* CSR_MPI_H_ */
