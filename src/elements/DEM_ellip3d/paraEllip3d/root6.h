///////////////////////////////////////////////////////////////////////////////////////////////////////
//                                   Code: ParaEllip3d-CFD                                           //
//                                 Author: Dr. Beichuan Yan                                          //
//                                  Email: beichuan.yan@colorado.edu                                 //
//                              Institute: University of Colorado Boulder                            //
///////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef ROOT6_H
#define ROOT6_H

#include "realtypes.h"
#include "Vec.h"

namespace dem {
  bool root6(REAL coef1[], REAL coef2[], Vec & v);
}
#endif
