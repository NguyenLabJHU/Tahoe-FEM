#ifndef CONST_H
#define CONST_H
#include "realtypes.h"
#include <fstream>

namespace dem { 

  // Pi
  extern const REAL Pi;

  // numerical EPS (NOT machine epsilon)
  extern const REAL EPS;

  // random number seed (NOT a constant)
  extern long idum;

  // output field width and precision
  extern const std::size_t OWID;
  extern const std::size_t OPREC;

  // other global variables
  extern std::ofstream debugInf;
  extern std::size_t iteration;
  extern REAL timeStep;

}
#endif
