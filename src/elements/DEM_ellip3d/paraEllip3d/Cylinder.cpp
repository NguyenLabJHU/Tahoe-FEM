#include "Cylinder.h"
#include "ran.h"
#include <iostream>

namespace dem {

void Cylinder::print() const{
    std::cout<<"radius="<<radius<<std::endl;
    std::cout<<"height="<<height<<std::endl;
    std::cout<<"center=";
    center.print();
}

Vec Cylinder::randomPoint() const{
    REAL temp1=ran(&idum);
    REAL temp2=ran(&idum);
    REAL temp3=ran(&idum);
    REAL z=(center.getz()+height/2)*temp1+(center.getz()-height/2)*(1-temp1);
    REAL theta=2*PI*temp2;
    REAL r=radius*temp3;
    REAL x=center.getx()+r*cos(theta);
    REAL y=center.gety()+r*sin(theta);
    return Vec(x,y,z); 
}

}
