// DEVELOPMENT
#ifndef _HEX_T_H_
#define _HEX_T_H_

#include "StringT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "CrystalLatticeT.h"

using namespace Tahoe;

class HEXT : public CrystalLatticeT 
{
public:
	HEXT(int nlsd,int nuca,dArrayT alat,
	     int which_rot,dArray2DT mat_rot,
	     double angle);

	~HEXT() { };

	HEXT(const HEXT& source);

	const dArrayT& GetVector_a(){ return a;};
	const dArrayT& GetVector_b(){ return b;};	
	const dArrayT& GetVector_c(){ return c;};

        const dArrayT& GetLatticeParameters();
        const dArray2DT& GetBasis();
	const dArray2DT& GetAxis();

};

#endif
