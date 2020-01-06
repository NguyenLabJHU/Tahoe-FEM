#ifndef CONST_H
#define CONST_H
#include <mpi.h>
#include "realtypes.h"
#include <cstddef>
#include <fstream>
#include <random>

namespace dem { 

  // Pi
  extern const REAL Pi;

  // numerical EPS (NOT machine epsilon)
  extern const REAL EPS;

  // algorithm tolerance
  extern const REAL TOL;

  // minimum particle number for computing granular stress and strain
  extern const int stressMinPtcl;

  // random number seed (NOT a constant)
  extern long idum;

  // declaration
  extern std::default_random_engine engine;
  extern std::uniform_real_distribution<double> ran11;

  // output field width and precision
  extern const std::size_t OWID;
  extern const std::size_t OPREC;

  // other global variables
  extern std::ofstream debugInf;
  extern std::size_t iteration;
  extern REAL timeStep;
  extern REAL timeAccrued;

#ifndef NDEBUG
  extern MPI_File overlapInf;
#endif

}
#endif
