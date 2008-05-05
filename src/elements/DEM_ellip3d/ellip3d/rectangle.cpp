#include <cmath>
#include <iostream>
#include "rectangle.h"
#include "ran.h"
using namespace std;

namespace dem {

extern long idum;

void rectangle::print() const{
    printf("%15.6Lf%15.6Lf%15.6Lf\n",width,length,height);
    printf("%15.6Lf%15.6Lf%15.6Lf\n",center.getx(),center.gety(),center.getz());
}

vec rectangle::randomPoint() const{
    long double temp1=ran(&idum);
    long double temp2=ran(&idum);
    long double temp3=ran(&idum);
    long double x=temp1*(center.getx()-width/2)+(1-temp1)*(center.getx()+width/2);
    long double y=temp2*(center.gety()-length/2)+(1-temp2)*(center.gety()+length/2);
    long double z=temp3*(center.getz()-height/2)+(1-temp3)*(center.getz()+height/2);
    return vec(x,y,z);
}

} // namespace dem ends
