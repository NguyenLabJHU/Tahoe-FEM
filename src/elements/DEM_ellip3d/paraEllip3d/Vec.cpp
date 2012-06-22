#include "Vec.h"
#include "parameter.h"
#include <iostream>

namespace dem {
  
  bool Vec::operator==(const Vec v){
    return x==v.x && y==v.y && z==v.z;
  }
  
  bool Vec::operator==(const REAL d){
    return x==d && y==d && z==d;
  }
  
  bool Vec::operator!=(const Vec v){
    return x!=v.x || y!=v.y || z!=v.z;
  }
  
  void Vec::operator+=(Vec v){
    x += v.x;
    y += v.y;
    z += v.z;
  }
  
  void Vec::operator-=(Vec v){
    x -= v.x;
    y -= v.y;
    z -= v.z;
  }
  
  void Vec::operator*=(REAL d){
    x *= d;
    y *= d;
    z *= d;
  }
  
  void Vec::operator/=(REAL d){
    x /= d;
    y /= d;
    z /= d;
  }
  
  Vec Vec::operator+(Vec v) const{		
    return Vec(x+v.x, y+v.y, z+v.z);
  }
  
  Vec Vec::operator-(Vec v) const{
    return Vec(x-v.x, y-v.y, z-v.z);
  }
  
  Vec Vec::operator*(Vec p) const{
    return Vec(y*p.z-z*p.y, z*p.x-x*p.z, x*p.y-y*p.x);
  }
  
  Vec Vec::operator*(REAL d) const{
    return Vec(x*d, y*d, z*d);
  }
  
  REAL Vec::operator%(Vec p) const{
    return (x*p.x + y*p.y + z*p.z);
  }
  
  void Vec::print() const{
    std::cout << "(" << x <<" "<< y << " " << z << ")" << std::endl;
  }
  
  Vec operator*(REAL d, Vec v){
    return Vec(v.getX()*d, v.getY()*d, v.getZ()*d);
  }

  Vec operator/(Vec v, REAL d){
    return Vec(v.getX()/d, v.getY()/d, v.getZ()/d);
  }
  
  REAL vfabs(Vec v){
    REAL x = v.getX();
    REAL y = v.getY();
    REAL z = v.getZ();
    return sqrt(x * x + y * y + z * z);
  }
  
  Vec vcos(Vec v){
    return Vec(cos(v.getX()), cos(v.getY()), cos(v.getZ()));
  }
  
  Vec vacos(Vec v){
    return Vec(acos(v.getX()), acos(v.getY()), acos(v.getZ()));
  }
  
  Vec operator-(Vec v){
    return -1.0*v;
  }
  
  Vec normalize(Vec v){
    REAL alf = vfabs(v);
    if (alf < EPS) // important, otherwise my cause numerical instability
      return v;
    return v/(vfabs(v));
  }
  
  Vec rotateVec(Vec v, Vec ang){
    REAL alf = vfabs(ang);
    if (alf < EPS) // important, otherwise my cause numerical instability
      return v;
    
    Vec nx = ang/alf;
    Vec vp = (v % nx) * nx;
    Vec vv = v - vp;
    
    REAL theta = atan(vfabs(vv) / vfabs(vp));
#ifndef NDEBUG
    debugInf<<"Vec.cpp: iter="<<iteration 
	      <<" alf="<<alf
	      <<" theta="<<theta<<std::endl;
#endif
    if (theta < EPS) // important, otherwise my cause numerical instability
      return v;    
    
    Vec ny=normalize(vv);
    Vec nz=normalize(nx*ny); // normalize, for higher precision
    REAL l=vfabs(vv);
    return l * sin(alf) * nz + l * cos(alf) * ny + vp;
  }
  
  REAL angle(Vec v1, Vec v2, Vec norm){
    //calculate the angle between v1 and v2 if rotating v1 in the plane
    //composed of v1 and v2 from itself to v2, the angle could be 0<alf<360
    //norm specify that the rotation must be around norm according to right hand rule,
    //even if the 180<alf<360
    REAL alf;
    Vec crs = v1 * v2;
    alf=asin(vfabs(crs) / vfabs(v1) / vfabs(v2) ); //0<alf<90;
    if(crs % norm > 0){ //0<=alf<=180
      if(v1 % v2 < 0)   //90<alf<180
	alf = PI-alf;
    }
    else{//180<alf<360
      if(v1 % v2 > 0)   //270<alf<360
	alf = 2 *PI - alf;
      else
	alf = PI + alf;
    }
    return alf;
  }
  
} // namespace dem
