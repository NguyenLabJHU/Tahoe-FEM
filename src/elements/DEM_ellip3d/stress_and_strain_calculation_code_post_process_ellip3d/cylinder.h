#ifndef CYLINDER_H
#define CYLINDER_H

#include "realtypes.h"
#include "parameter.h"
#include "vec.h"
#include "shape.h"

namespace dem {

class cylinder:public shape{
 public:
    cylinder()
	:radius(0),height(0),center(0)
	{}

    cylinder(REAL r, REAL h, vec c)
	:radius(r),height(h),center(c)
	{}

    cylinder(const cylinder &cy){
	radius=cy.radius;
	height=cy.height;
	center=cy.center;
    }
    
    REAL getRadius() const {return radius;}
    REAL getHeight() const {return height;}
    REAL getVolume() const {return PI*radius*radius*height;}
    vec  getCenter() const {return center;}

    void setRadius(REAL r) {radius = r;}
    void setHeight(REAL h) {height = h;}
    void setCenter(vec v) {center=v;}

    vec  randomPoint() const;
    void print() const;
    
 private:
    REAL radius;
    REAL height;
    vec  center;
};

}

#endif
