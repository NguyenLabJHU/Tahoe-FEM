#include "const.h"

namespace dem { 

  // Pi
  const REAL Pi   = 3.141592653589;

  // numerical EPS (NOT machine epsilon)
  const REAL EPS  = 1.0E-12;

  // minimum particle number for computing granular stress and strain
  // mininum 5 to construct initial simplex
  const int stressMinPtcl = 5;

  // algorithm tolerance
  const REAL TOL  = 1.0E-6;

  // random number seed (Not a constant)
  long idum       = -1;

  // definition
  std::default_random_engine engine;
  std::uniform_real_distribution<double> ran11(0.0, 1.0);

  // output field width and precision
  const std::size_t OWID  = 15;   // output width
  const std::size_t OPREC = 6;    // output precision, number of digits after decimal dot

  // other global variables
  std::ofstream debugInf;         // debug info, only root process prints to debugInf
  std::size_t iteration;          // iteration number
  REAL timeStep;                  // time step
  REAL timeAccrued;               // accurued time
#ifndef NDEBUG
  MPI_File overlapInf;            // contact overlap info, parallel IO
#endif

}
