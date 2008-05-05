#ifndef SHAPE_H
#define SHAPE_H
#include "vec.h"

namespace dem {

class shape{ // abstract class
public:
    virtual vec get_center() const =0;
    virtual long double get_volume() const =0;
    virtual vec randomPoint() const =0;
    virtual void  print() const =0;
    virtual        ~shape() {}; // base class needs a virtual destructor.
};

} // namespace dem ends

#endif
