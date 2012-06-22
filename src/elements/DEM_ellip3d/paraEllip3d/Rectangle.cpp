#include "Rectangle.h"
#include "parameter.h"
#include "ran.h"
#include <iostream>

namespace dem {

  void Rectangle::print() const{
    std::cout << v1.getX() << ' ' << v1.getY() << ' ' << v1.getZ() << ' '
	      << v2.getX() << ' ' << v2.getY() << ' ' << v2.getZ() << std::endl;
  }
  
  Vec Rectangle::randomPoint() const{
    REAL tmp1 = ran(&idum);
    REAL tmp2 = ran(&idum);
    REAL tmp3 = ran(&idum);
    REAL x = tmp1*(center.getX() - dimx/2) + (1-tmp1)*(center.getX() + dimx/2);
    REAL y = tmp2*(center.getY() - dimy/2) + (1-tmp2)*(center.getY() + dimy/2);
    REAL z = tmp3*(center.getZ() - dimz/2) + (1-tmp3)*(center.getZ() + dimz/2);
    return Vec(x,y,z);
  }
  
} // namespace dem
