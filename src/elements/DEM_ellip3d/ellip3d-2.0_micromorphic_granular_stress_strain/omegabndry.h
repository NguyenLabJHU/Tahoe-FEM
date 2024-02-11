#ifndef OMEGABNDRY_H
#define OMEGABNDRY_H

#include "realtypes.h"
#include "matrix.h"
#include "microrve.h"


namespace dem {


// this class is for the boundary faces of Omega which is the domain to calculate micromophic stress
class omegabndry{

friend class assembly;

private:
    microrve * RVE;	// rve is the RVE to which this boundary face belongs
    matrix n_bndry;	// unit normal vector on boundary of this RVE, column vector
    REAL area;		// area of this boundary face

    int num;		// 1 means x boundary, 2 means y boundary, 3 means z boundary


public:
    omegabndry(microrve* rve, matrix n, REAL a, int num_t){RVE = rve; n_bndry = n; area = a; num=num_t;}	

}; // class omegabndry



} // namespace dem



#endif
