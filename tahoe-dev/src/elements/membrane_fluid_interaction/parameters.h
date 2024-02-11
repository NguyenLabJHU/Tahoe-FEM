#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "realtypes.h"

namespace memFluid{

    // these material properties in matlab files are lower case
    extern const REAL MU_FLUID;
    extern const REAL RHO;

    extern const REAL KC;
    extern const REAL GAMMA;
    extern const REAL NU;
    extern const REAL LAMBDA_MEMBRANE;
    extern const REAL MU_MEMBRANE;

    extern const REAL MU_F_MINUS;
    extern const REAL MU_F_PLUS;

    extern REAL L_PLUS;
    extern REAL L_MINUS;

    extern REAL MU_PLUS;
    extern REAL MU_MINUS;
////////  materials end

    // floor position
    extern REAL YP;
    extern const int REPULSION1;
    extern const int REPULSION2;

    extern const REAL PI;

    // output field width and precision
    extern const int OWID;
    extern const int OPREC;

    extern const REAL sparseTHRESHOLD;


} // end of memFluid

#endif
