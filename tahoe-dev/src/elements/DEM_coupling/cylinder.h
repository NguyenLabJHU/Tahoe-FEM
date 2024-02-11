#ifndef CYLINDER_H
#define CYLINDER_H

#include "const.h"
#include "vec.h"
#include "shape.h"

namespace dem {

class cylinder:public shape{
 public:
    cylinder()
	:radius(0),height(0),center(0)
	{}

    cylinder(long double r, long double h, vec c)
	:radius(r),height(h),center(c)
	{}

    cylinder(const cylinder &cy){
	radius=cy.radius;
	height=cy.height;
	center=cy.center;
    }
    
    long double get_radius() const {return radius;}
    long double get_height() const {return height;}
    long double get_volume() const {return PI*radius*radius*height;}
    vec  get_center() const {return center;}

    void set_radius(long double r) {radius = r;}
    void set_height(long double h) {height = h;}
    void set_center(vec v) {center=v;}

    vec  randomPoint() const;
    void print() const;
    
 private:
    long double radius;
    long double height;
    vec  center;
};

}

#endif
