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
    
    REAL get_radius() const {return radius;}
    REAL get_height() const {return height;}
    REAL get_volume() const {return PI*radius*radius*height;}
    vec  get_center() const {return center;}

    void set_radius(REAL r) {radius = r;}
    void set_height(REAL h) {height = h;}
    void set_center(vec v) {center=v;}

    vec  randomPoint() const;
    void print() const;
    
 private:
    REAL radius;
    REAL height;
    vec  center;
};

}

#endif
