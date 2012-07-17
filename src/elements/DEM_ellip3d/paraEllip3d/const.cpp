#include "const.h"

namespace dem { 

// Pi
const REAL Pi   = 3.141592653589;

// numerical EPS (NOT machine epsilon)
const REAL EPS  = 1.0E-12;

// random number seed (Not a constant)
long idum       = -1;

// output field width and precision
const int OWID  = 16;   // 20, output width
const int OPREC = 6;    // 10, output precision, number of digits after decimal dot

// other global variables
std::ofstream debugInf; // debug info
int iteration;          // iteration number

}
