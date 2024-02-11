#include "rectangle.h"
#include "parameter.h"
#include "ran.h"
#include <iostream>

namespace dem {

  void rectangle::print() const{
    std::cout << "dimensions: " << dimx << " " << dimy << " " << dimz << std::endl;
    std::cout << "center: " << center.getx() << " " << center.gety() << " " << center.getz() << std::endl;
    std::cout << "reference: " << v1.getx() << " " << v1.gety() << " " << v1.getz() << std::endl;
  }
  
  vec rectangle::randomPoint() const{
    REAL tmp1=ran(&idum);
    REAL tmp2=ran(&idum);
    REAL tmp3=ran(&idum);
    REAL x=tmp1*(center.getx()-dimx/2)+(1-tmp1)*(center.getx()+dimx/2);
    REAL y=tmp2*(center.gety()-dimy/2)+(1-tmp2)*(center.gety()+dimy/2);
    REAL z=tmp3*(center.getz()-dimz/2)+(1-tmp3)*(center.getz()+dimz/2);
    return vec(x,y,z);
  }
  
} // namespace dem
