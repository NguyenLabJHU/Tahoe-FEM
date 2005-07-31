#ifndef _SPOOLES_MPI_H_
#define _SPOOLES_MPI_H_

#include "SPOOLES.h"
#include "mpi.h"

#ifdef __cplusplus
extern "C" {
#endif
/* MPI driver provided by in SPOOLES documentation */
extern int LU_MPI_driver(int msglvl, const char* message_file, int matrix_type,
	int symmetry_flag, int pivoting_flag, int seed, int num_eq, int num_row,
	int* r_rhs, double* rhs2out, int num_entries, int* r, int* c, double* v, MPI_Comm* comm);
#ifdef __cplusplus
}
#endif
#endif /* _SPOOLES_MPI_H_ */
