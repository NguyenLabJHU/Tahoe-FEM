///////////////////////////////////////////////////////////////////////////////////////////////////////
//                                   Code: ParaEllip3d-CFD                                           //
//                                 Author: Dr. Beichuan Yan                                          //
//                                  Email: beichuan.yan@colorado.edu                                 //
//                              Institute: University of Colorado Boulder                            //
///////////////////////////////////////////////////////////////////////////////////////////////////////

#include "Cylinder.h"
#include "ran.h"
#include <iostream>

namespace dem {

  void Cylinder::print() const{
    debugInf << "radius=" << radius << std::endl;
    debugInf << "height=" << height << std::endl;
    debugInf << "center=";
    center.print(debugInf);
  }

  Vec Cylinder::randomPoint() const{
    REAL rand1 = ran11(engine);
    REAL rand2 = ran11(engine);
    REAL rand3 = ran11(engine);
    REAL z = (center.getZ() + height/2)*rand1 + (center.getZ() - height/2)*(1 - rand1);
    REAL theta = 2*Pi*rand2;
    REAL r = radius*rand3;
    REAL x = center.getX() + r*cos(theta);
    REAL y = center.getY() + r*sin(theta);
    return Vec(x,y,z); 
  }

}
