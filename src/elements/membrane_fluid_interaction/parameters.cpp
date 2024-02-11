#include "parameters.h"

namespace memFluid{


    // these material properties in matlab files are lower case
    const REAL MU_FLUID = 1e-10;
    const REAL RHO = 2e-7;

    const REAL KC = 0;
    const REAL GAMMA = 1e-6;
    const REAL NU = 0;
    const REAL LAMBDA_MEMBRANE = 0;
    const REAL MU_MEMBRANE = 0;

    const REAL MU_F_MINUS = MU_FLUID;
    const REAL MU_F_PLUS = MU_FLUID;

    REAL L_PLUS;
    REAL L_MINUS;

    REAL MU_PLUS;
    REAL MU_MINUS;
////////  materials end

    // floor position
    REAL YP;
    const int REPULSION1 = 100;
    const int REPULSION2 = 20;

    const REAL PI = 3.141592653589;

    // output width and precision
    const int OWID        = 16;      // 20, output width
    const int OPREC       = 6;       // 10, output precision, number of digits after decimal dot

    const REAL sparseTHRESHOLD = 0.0000000001;


} // end of memFluid
