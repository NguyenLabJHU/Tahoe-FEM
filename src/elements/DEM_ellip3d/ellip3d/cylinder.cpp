#include "cylinder.h"
#include "ran.h"
#include <iostream>
#include <cmath>

using namespace std;

namespace dem {

void cylinder::print() const{
    cout<<"radius="<<radius<<endl;
    cout<<"height="<<height<<endl;
    cout<<"center=";
    center.print();
}

vec cylinder::randomPoint() const{
    REAL temp1=ran(&idum);
    REAL temp2=ran(&idum);
    REAL temp3=ran(&idum);
    REAL z=(center.getz()+height/2)*temp1+(center.getz()-height/2)*(1-temp1);
    REAL theta=2*PI*temp2;
    REAL r=radius*temp3;
    REAL x=center.getx()+r*cosl(theta);
    REAL y=center.gety()+r*sinl(theta);
    return vec(x,y,z); 
}

}
