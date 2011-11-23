#ifndef GRADATION_H
#define GRADATION_H

#include "realtypes.h"
#include <vector>

namespace dem {

class gradation{
public:
    int rorc;
    REAL dimn;
    REAL ratio_ba;
    REAL ratio_ca;
    int seivenum;
    std::vector<REAL> percent;
    std::vector<REAL> ptclsize;

    gradation()
	:rorc(),dimn(),ratio_ba(),ratio_ca(),seivenum(),percent(),ptclsize()
	{};

    gradation(int _rc, REAL _di, REAL _ba, REAL _ca,
	      int _sn, std::vector<REAL>& _v1, std::vector<REAL>& _v2)
	:rorc(_rc), dimn(_di), ratio_ba(_ba), ratio_ca(_ca),
	seivenum(_sn), percent(_v1), ptclsize(_v2)
	{};

};

}

#endif
