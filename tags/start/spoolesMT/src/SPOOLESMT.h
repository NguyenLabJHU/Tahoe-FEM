#ifndef _SPOOLES_MT_H_
#define _SPOOLES_MT_H_

#include "SPOOLES.h"

#ifdef __cplusplus
extern "C" {
#endif
/* MT driver provided by in SPOOLES documentation */
extern int  LU_MT_driver(int msg_lvl, const char* message_file, 
	    int matrix_type, int symmetry_flag, int pivoting_flag, 
            int rand_seed, int num_eq, double* rhs2out, int num_entries, 
			 int* r, int* c, double* v, int nthread);

#ifdef __cplusplus
}
#endif
#endif /* _SPOOLES_MT_H_ */





