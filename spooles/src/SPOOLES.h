#ifndef _SPOOLES__H_
#define _SPOOLES__H_

#define SPOOLES_INDICES_ONLY    0
#define SPOOLES_REAL            1
#define SPOOLES_COMPLEX         2

#define SPOOLES_SYMMETRIC    0
#define SPOOLES_HERMITIAN    1
#define SPOOLES_NONSYMMETRIC 2

#define SPOOLES_NO_PIVOTING  0
#define SPOOLES_PIVOTING     1

#define SPOOLES_BY_ROWS      1
#define SPOOLES_BY_COLUMNS   2

#ifdef __cplusplus
extern "C" {
#endif

/* serial driver provided by in SPOOLES documentation */
extern int LU_serial_driver(int msglvl, const char* message_file, int matrix_type,
	int symmetry_flag, int pivoting_flag, int seed, int num_eq, double* rhs2out, 
	int num_entries, int* r, int* c, double* v);

#ifdef __cplusplus
}
#endif
#endif
