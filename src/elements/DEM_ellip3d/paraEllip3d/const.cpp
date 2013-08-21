#include "const.h"

namespace dem { 

// Pi
const REAL Pi   = 3.141592653589;

// numerical EPS (NOT machine epsilon)
const REAL EPS  = 1.0E-12;

// random number seed (Not a constant)
long idum       = -1;

// output field width and precision
const std::size_t OWID  = 16;   // 20, output width
const std::size_t OPREC = 6;    // 10, output precision, number of digits after decimal dot

// other global variables
std::ofstream debugInf;         // debug info
std::size_t iteration;          // iteration number

}
