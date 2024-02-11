#include <iostream>
#include <cmath>
#include "ran.h"
#include "cylinder.h"
using namespace std;

namespace dem {

extern long idum;

void cylinder::print() const{
    cout<<"radius="<<radius<<endl;
    cout<<"height="<<height<<endl;
    cout<<"center=";
    center.print();
}

vec cylinder::randomPoint() const{
    long double temp1=ran(&idum);
    long double temp2=ran(&idum);
    long double temp3=ran(&idum);
    long double z=(center.getz()+height/2)*temp1+(center.getz()-height/2)*(1-temp1);
    long double theta=2*PI*temp2;
    long double r=radius*temp3;
    long double x=center.getx()+r*cosl(theta);
    long double y=center.gety()+r*sinl(theta);
    return vec(x,y,z); 
}

}
