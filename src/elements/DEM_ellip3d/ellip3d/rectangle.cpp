#include "rectangle.h"
#include "parameter.h"
#include "ran.h"
#include <iostream>
#include <cmath>
#include <cstdio>
using namespace std;

namespace dem {

void rectangle::print() const{
    printf("%15.6Lf%15.6Lf%15.6Lf\n",width,length,height);
    printf("%15.6Lf%15.6Lf%15.6Lf\n",center.getx(),center.gety(),center.getz());
}

vec rectangle::randomPoint() const{
    REAL temp1=ran(&idum);
    REAL temp2=ran(&idum);
    REAL temp3=ran(&idum);
    REAL x=temp1*(center.getx()-width/2)+(1-temp1)*(center.getx()+width/2);
    REAL y=temp2*(center.gety()-length/2)+(1-temp2)*(center.gety()+length/2);
    REAL z=temp3*(center.getz()-height/2)+(1-temp3)*(center.getz()+height/2);
    return vec(x,y,z);
}

} // namespace dem ends
