#ifndef VEC_H
#define VEC_H

#include "realtypes.h"

namespace dem {

class vec{
 public:
    vec();
    vec(REAL d);
    vec(REAL _x, REAL _y, REAL _z);
    
    REAL getx() const;
    REAL gety() const;
    REAL getz() const;
    void setx(REAL _x);
    void sety(REAL _y);
    void setz(REAL _z);
    
    bool operator==(const vec v){
	return x==v.x && y==v.y && z==v.z;
    }
    
    bool operator==(const REAL d){
	return x==d && y==d && z==d;
    }
    
    bool operator!=(const vec v){
	return x!=v.x || y!=v.y || z!=v.z;
    }
    
    void operator+=(const vec v);
    void operator-=(const vec v);
    void operator*=(REAL d);
    void operator/=(REAL d);
    vec operator+(vec v) const;
    vec operator-(vec v) const;
    vec operator*(vec p) const;   // find the cross product of this vector and p
    vec operator*(REAL d) const;
    REAL operator%(vec p) const;// find the dot product of this and p
    void print() const;

 private:
    REAL x;
    REAL y;
    REAL z;
};

// Non-member functions
vec operator*(REAL d, vec v);
vec operator/(vec v, REAL d);
vec operator-(vec v);
REAL vfabs(vec v);
vec vcos(vec v);
vec vacos(vec v);
vec rotateVec(vec v, vec alf);    // find the exact vector after v is rotated alf in space
vec normalize(vec v);
/*calculate the angle between v1 and v2 if rotating v1 in the plane
composed of v1 and v2 from itself to v2, the angle could be 0<alf<360
norm specify that the rotation must be around norm according to right hand rule,
even if 180<alf<360
*/
REAL angle(vec v1, vec v2, vec norm); 

} // namespace dem ends

#endif
