#ifndef RECTANGLE_H
#define RECTANGLE_H

#include "realtypes.h"
#include "vec.h"
#include "shape.h"

namespace dem {

class rectangle:public shape{
 public:
    rectangle()
	:width(0),length(0),height(0),center(0)
	{}

    rectangle(REAL w, REAL l, REAL h, vec c)
	:width(w),length(l),height(h),center(c)
	{}

    rectangle(const rectangle &rec) {
	width=rec.width;
	length=rec.length;
	height=rec.height;
	center=rec.center;
    }
    
    REAL get_width() const {return width;}
    REAL get_length() const {return length;}
    REAL get_height() const {return height;}
    REAL get_volume() const {return width*length*height;}
    vec  get_center() const {return center;}

    void set_width(REAL w) {width=w;}
    void set_length(REAL l) {length=l;}
    void set_height(REAL h) {height=h;}
    void set_center(vec v) {center=v;}

    vec  randomPoint() const;
    void print() const;
    
 private:
    REAL width;
    REAL length;
    REAL height;
    vec  center;
};

} // namespace dem ends

#endif
