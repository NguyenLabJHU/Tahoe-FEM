#ifndef OMEGABNDRY_H
#define OMEGABNDRY_H

#include "realtypes.h"
#include "matrix.h"
#include "microrve"


namespace dem {


// this class is for the boundary faces of Omega which is the domain to calculate micromophic stress
class omegabndry{

friend class assembly;

private:
    microrve * RVE;	// rve is the RVE to which this boundary face belongs
    matrix n_bndry;	// unit normal vector on boundary of this RVE, column vector
    REAL area;		// area of this boundary face


public:
    omegabndry(microrve* rve, matrix n, REAL d){RVE = rve; n_bndry = n; area = d*d;}	

}; // class omegabndry



} // namespace dem



#endif
