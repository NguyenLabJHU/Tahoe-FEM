#ifndef RECTANGLE_H
#define RECTANGLE_H

#include "vec.h"
#include "shape.h"

namespace dem {

class rectangle:public shape{
 public:
    rectangle()
	:width(0),length(0),height(0),center(0)
	{}

    rectangle(long double w, long double l, long double h, vec c)
	:width(w),length(l),height(h),center(c)
	{}

    rectangle(const rectangle &rec) {
	width=rec.width;
	length=rec.length;
	height=rec.height;
	center=rec.center;
    }
    
    long double get_width() const {return width;}
    long double get_length() const {return length;}
    long double get_height() const {return height;}
    long double get_volume() const {return width*length*height;}
    vec  get_center() const {return center;}

    void set_width(long double w) {width=w;}
    void set_length(long double l) {length=l;}
    void set_height(long double h) {height=h;}
    void set_center(vec v) {center=v;}

    vec  randomPoint() const;
    void print() const;
    
 private:
    long double width;
    long double length;
    long double height;
    vec  center;
};

} // namespace dem ends

#endif
