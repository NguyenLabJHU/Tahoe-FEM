#ifndef ARRAY_H
#define ARRAY_H

#include "realtypes.h"

namespace memFluid{


class array {

private:
    REAL ptx;
    REAL pty;

public:
    array(REAL val1, REAL val2): ptx(val1), pty(val2) {}
    REAL getX() const {return ptx;}
    REAL getY() const {return pty;}

    // non member functions
    friend /*inline */bool operator < (const array&, const array&);

};



} // end memFluid


#endif
