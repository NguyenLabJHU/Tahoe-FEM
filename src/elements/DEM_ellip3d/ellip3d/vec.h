#ifndef VEC_H
#define VEC_H

namespace dem {

class vec{
 public:
    vec();
    vec(long double d);
    vec(long double _x, long double _y, long double _z);
    
    long double getx() const;
    long double gety() const;
    long double getz() const;
    void setx(long double _x);
    void sety(long double _y);
    void setz(long double _z);
    
    bool operator==(const vec v){
	return x==v.x && y==v.y && z==v.z;
    }
    
    bool operator==(const long double d){
	return x==d && y==d && z==d;
    }
    
    bool operator!=(const vec v){
	return x!=v.x || y!=v.y || z!=v.z;
    }
    
    void operator+=(const vec v);
    void operator-=(const vec v);
    void operator*=(long double d);
    void operator/=(long double d);
    vec operator+(vec v) const;
    vec operator-(vec v) const;
    vec operator*(vec p) const;   // find the cross product of this vector and p
    vec operator*(long double d) const;
    long double operator%(vec p) const;// find the dot product of this and p
    void print() const;
    
 private:
    long double x;
    long double y;
    long double z;
};

// Non-member functions
vec operator*(long double d, vec v);
vec operator/(vec v, long double d);
vec operator-(vec v);
long double vfabsl(vec v);
vec vcosl(vec v);
vec vacosl(vec v);
vec rotateVec(vec v, vec alf);    // find the exact vector after v is rotated alf in space
vec normalize(vec v);
/*calculate the angle between v1 and v2 if rotating v1 in the plane
composed of v1 and v2 from itself to v2, the angle could be 0<alf<360
norm specify that the rotation must be around norm according to right hand rule,
even if 180<alf<360
*/
long double angle(vec v1, vec v2, vec norm); 

} // namespace dem ends

#endif
