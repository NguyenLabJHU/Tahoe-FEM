#include "Cylinder.h"
#include "ran.h"
#include <iostream>

namespace dem {

void Cylinder::print() const{
    std::cout << "radius=" << radius << std::endl;
    std::cout << "height=" << height << std::endl;
    std::cout << "center=";
    center.print();
}

Vec Cylinder::randomPoint() const{
    REAL tmp1 = ran(&idum);
    REAL tmp2 = ran(&idum);
    REAL tmp3 = ran(&idum);
    REAL z = (center.getZ() + height/2)*tmp1 + (center.getZ() - height/2)*(1 - tmp1);
    REAL theta = 2*PI*tmp2;
    REAL r = radius*tmp3;
    REAL x = center.getX() + r*cos(theta);
    REAL y = center.getY() + r*sin(theta);
    return Vec(x,y,z); 
}

}
