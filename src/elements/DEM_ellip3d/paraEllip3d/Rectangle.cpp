#include "Rectangle.h"
#include "const.h"
#include "ran.h"
#include <iostream>

namespace dem {

  void Rectangle::print() const {
    debugInf << v1.getX() << ' ' << v1.getY() << ' ' << v1.getZ() << ' '
	     << v2.getX() << ' ' << v2.getY() << ' ' << v2.getZ() << std::endl;
  }
  
  Vec Rectangle::randomPoint() const {
    REAL rand1 = ran(&idum);
    REAL rand2 = ran(&idum);
    REAL rand3 = ran(&idum);
    REAL x = rand1*(center.getX() - dimX/2) + (1-rand1)*(center.getX() + dimX/2);
    REAL y = rand2*(center.getY() - dimY/2) + (1-rand2)*(center.getY() + dimY/2);
    REAL z = rand3*(center.getZ() - dimZ/2) + (1-rand3)*(center.getZ() + dimZ/2);
    return Vec(x,y,z);
  }
  
  // the order 1, 2, 3, 4, 5, 6 corresponds to surface at -x, +x, -y, +y, -z, +z, respectively.
  Vec Rectangle::getSurfaceCenter(int i) const {
    Vec sCenter;
    switch (i) {
    case 1:
      sCenter = Vec(center.getX() - dimX/2, center.getY(), center.getZ());
      break;
    case 2:
      sCenter = Vec(center.getX() + dimX/2, center.getY(), center.getZ());
      break;
    case 3:
      sCenter = Vec(center.getX(), center.getY() - dimY/2, center.getZ());
      break;
    case 4:
      sCenter = Vec(center.getX(), center.getY() + dimY/2, center.getZ());
      break;
    case 5:
      sCenter = Vec(center.getX(), center.getY(), center.getZ() - dimZ/2);
      break;
    case 6:
      sCenter = Vec(center.getX(), center.getY(), center.getZ() + dimZ/2);
      break;
    }
    return sCenter;
  }

  // the vertices follow FEM convention
  Vec Rectangle::getVertice(int i) const {
    Vec vertice;
    switch (i) {
    case 1:
      vertice = v1;
      break;
    case 2:
      vertice = Vec(v1.getX() + dimX, v1.getY(), v1.getZ());
      break;
    case 3:
      vertice = Vec(v1.getX() + dimX, v1.getY() + dimY, v1.getZ());
      break;
    case 4:
      vertice = Vec(v1.getX(), v1.getY() + dimY, v1.getZ());
      break;
    case 5:
      vertice = Vec(v2.getX() - dimX, v2.getY() - dimY, v2.getZ());
      break;
    case 6:
      vertice = Vec(v2.getX(), v2.getY() - dimY, v2.getZ());
      break;
    case 7:
      vertice = v2;
      break;
    case 8:
      vertice = Vec(v2.getX() - dimX, v2.getY(), v2.getZ());
      break;
    }
    return vertice;
  }

} // namespace dem
