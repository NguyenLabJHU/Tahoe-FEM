#include "vec.h"
#include "parameter.h"
#include <iostream>

namespace dem {
  
  bool vec::operator==(const vec v){
    return x==v.x && y==v.y && z==v.z;
  }
  
  bool vec::operator==(const REAL d){
    return x==d && y==d && z==d;
  }
  
  bool vec::operator!=(const vec v){
    return x!=v.x || y!=v.y || z!=v.z;
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
  
  void vec::operator*=(REAL d){
    x *= d;
    y *= d;
    z *= d;
  }
  
  void vec::operator/=(REAL d){
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
  
  vec vec::operator*(REAL d) const{
    return vec(x*d, y*d, z*d);
  }
  
  REAL vec::operator%(vec p) const{
    return (x*p.x + y*p.y + z*p.z);
  }
  
  void vec::print() const{
    std::cout << "(" << x <<" "<< y << " " << z << ")" << std::endl;
  }
  
  vec operator*(REAL d, vec v){
    return vec(v.getx()*d, v.gety()*d, v.getz()*d);
  }

  vec operator/(vec v, REAL d){
    return vec(v.getx()/d, v.gety()/d, v.getz()/d);
  }
  
  REAL vfabs(vec v){
    REAL x = v.getx();
    REAL y = v.gety();
    REAL z = v.getz();
    return sqrt(x * x + y * y + z * z);
  }
  
  vec vcos(vec v){
    return vec(cos(v.getx()), cos(v.gety()), cos(v.getz()));
  }
  
  vec vacos(vec v){
    return vec(acos(v.getx()), acos(v.gety()), acos(v.getz()));
  }
  
  vec operator-(vec v){
    return -1.0*v;
  }
  
  vec normalize(vec v){
    REAL alf = vfabs(v);
    if (alf < EPS) // important, otherwise my cause numerical instability
      return v;
    return v/(vfabs(v));
  }
  
  vec rotateVec(vec v, vec ang){
    REAL alf = vfabs(ang);
    if (alf < EPS) // important, otherwise my cause numerical instability
      return v;
    
    vec nx = ang/alf;
    vec vp = (v % nx) * nx;
    vec vv = v - vp;
    
    REAL theta = atan(vfabs(vv) / vfabs(vp));
#ifndef NDEBUG
    g_debuginf<<"vec.cpp: g_iteration="<<g_iteration 
	      <<" alf="<<alf
	      <<" theta="<<theta<<std::endl;
#endif
    if (theta < EPS) // important, otherwise my cause numerical instability
      return v;    
    
    vec ny=normalize(vv);
    vec nz=normalize(nx*ny); // normalize, for higher precision
    REAL l=vfabs(vv);
    return l * sin(alf) * nz + l * cos(alf) * ny + vp;
  }
  
  REAL angle(vec v1, vec v2, vec norm){
    //calculate the angle between v1 and v2 if rotating v1 in the plane
    //composed of v1 and v2 from itself to v2, the angle could be 0<alf<360
    //norm specify that the rotation must be around norm according to right hand rule,
    //even if the 180<alf<360
    REAL alf;
    vec crs = v1 * v2;
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
