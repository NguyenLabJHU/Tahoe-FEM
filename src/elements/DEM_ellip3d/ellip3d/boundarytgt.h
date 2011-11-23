#ifndef BOUNDARYTGT_H
#define BOUNDARYTGT_H

#include "realtypes.h"
#include "vec.h"

namespace dem {

class boundarytgt{
public:
    int  ptcl;
    vec  TgtForce;
    vec  TgtDisp;
    bool TgtLoading;
    vec  TgtDispStart;
    REAL TgtPeak;

    boundarytgt();
    boundarytgt(int _ptcl, vec _v1, vec _v2, bool _b, vec _v3, REAL _tp)
	:ptcl(_ptcl), 
	 TgtForce(_v1), 
	 TgtDisp(_v2), 
	 TgtLoading(_b),
	 TgtDispStart(_v3),
	 TgtPeak(_tp)
	{};
};

}

#endif
