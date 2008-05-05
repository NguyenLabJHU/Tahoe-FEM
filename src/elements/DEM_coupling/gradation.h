#ifndef GRADATION_H
#define GRADATION_H
#include <vector>

namespace dem {

class gradation{
public:
    int rorc;
    long double dimn;
    long double ratio_ba;
    long double ratio_ca;
    int seivenum;
    std::vector<long double> percent;
    std::vector<long double> ptclsize;

    gradation()
	:rorc(),dimn(),ratio_ba(),ratio_ca(),seivenum(),percent(),ptclsize()
	{};

    gradation(int _rc, long double _di, long double _ba, long double _ca,
	      int _sn, std::vector<long double>& _v1, std::vector<long double>& _v2)
	:rorc(_rc), dimn(_di), ratio_ba(_ba), ratio_ca(_ca),
	seivenum(_sn), percent(_v1), ptclsize(_v2)
	{};

};

}

#endif
