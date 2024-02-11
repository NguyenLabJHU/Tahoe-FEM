#include "Vec.h"
#include <iostream>
  
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
  
  Vec Vec::operator%(Vec p) const{
    return Vec(y*p.z-z*p.y, z*p.x-x*p.z, x*p.y-y*p.x);
  }
  
  Vec Vec::operator*(REAL d) const{
    return Vec(x*d, y*d, z*d);
  }
  
  REAL Vec::operator*(Vec p) const{
    return (x*p.x + y*p.y + z*p.z);
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
