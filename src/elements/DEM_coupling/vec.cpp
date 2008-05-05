#include <cmath>
#include <iostream>
#include "vec.h"
#include "const.h"

using namespace std;

namespace dem {

vec::vec(){
    x=0;
    y=0;
    z=0;
}


vec::vec(long double d){
    x=d;
    y=d;
    z=d;
}


vec::vec(long double _x, long double _y, long double _z){
    x=_x;
    y=_y;
    z=_z;
}


long double vec::getx() const{
    return x;
}


long double vec::gety() const{
    return y;
}


long double vec::getz() const{
    return z;
}


void vec::setx(long double _x){
    x=_x;
}


void vec::sety(long double _y){
    y=_y;
}


void vec::setz(long double _z){
    z=_z;
}


void vec::operator+=(vec v){
    x += v.x;
    y += v.y;
    z += v.z;
}


void vec::operator-=(vec v){
    x -= v.x;
    y -= v.y;
    z -= v.z;
}


void vec::operator*=(long double d){
    x *= d;
    y *= d;
    z *= d;
}


void vec::operator/=(long double d){
    x /= d;
    y /= d;
    z /= d;
}


vec vec::operator+(vec v) const{		
    return vec(x+v.x, y+v.y, z+v.z);
}


vec vec::operator-(vec v) const{
    return vec(x-v.x, y-v.y, z-v.z);
}


vec vec::operator*(vec p) const{
    return vec(y*p.z-z*p.y, z*p.x-x*p.z, x*p.y-y*p.x);
}


vec vec::operator*(long double d) const{
    return vec(x*d, y*d, z*d);
}


long double vec::operator%(vec p) const{
    return (x*p.x + y*p.y + z*p.z);
}


void vec::print() const{
	cout<<'('<<x<<' '<<y<<' '<<z<<')'<<endl;
}


vec operator*(long double d, vec v){
    return vec(v.getx()*d, v.gety()*d, v.getz()*d);
}


vec operator/(vec v, long double d){
    return vec(v.getx()/d, v.gety()/d, v.getz()/d);
}


long double vfabsl(vec v){
    return sqrtl(v.getx()*v.getx()+v.gety()*v.gety()+v.getz()*v.getz());
}


vec vcosl(vec v){
    return vec(cosl(v.getx()), cosl(v.gety()), cosl(v.getz()));
}


vec vacosl(vec v){
    return vec(acosl(v.getx()), acosl(v.gety()), acosl(v.getz()));
}


vec operator-(vec v){
	return -1*v;
}


vec normalize(vec v){
    return v/(vfabsl(v));
}


vec rotateVec(vec v, vec ang){
    long double alf=vfabsl(ang);
    if (alf<PREC) 
	return v;

    vec nx=ang/alf;
    vec vp=(v%nx)*nx;
    vec vv=v-vp;

    long double theta=atanl(vfabsl(vv)/vfabsl(vp));
    if (theta<PREC) 
	return v;    

    vec ny=normalize(vv);
    vec nz=normalize(nx*ny); // normalize, for better precision
    long double l=vfabsl(vv);
    return l*sinl(alf)*nz + l*cosl(alf)*ny + vp;
}


long double angle(vec v1, vec v2, vec norm){
    //calculate the angle between v1 and v2 if rotating v1 in the plane
    //composed of v1 and v2 from itself to v2, the angle could be 0<alf<360
    //norm specify that the rotation must be around norm according to right hand rule,
    //even if the 180<alf<360
    long double alf;
    vec crs=v1*v2;
    alf=asinl(vfabsl(crs)/vfabsl(v1)/vfabsl(v2));//0<alf<90;
    if(crs%norm>0){//0<=alf<=180
	if(v1%v2<0)//90<alf<180
	    alf=PI-alf;
    }
    else{//180<alf<360
	if(v1%v2>0)//270<alf<360
	    alf=2*PI-alf;
	else
	    alf=PI+alf;
    }
    return alf;
}

} // namespace dem ends
