// DEVELOPMENT
#ifndef _DIA_T_H_
#define _DIA_T_H_

#include "StringT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "CrystalLatticeT.h"

using namespace Tahoe;

class DIAT : public CrystalLatticeT 
{
public:
	DIAT(int nlsd,int nuca,dArrayT alat,
	     int which_rot,dArray2DT mat_rot,
	     double angle);

	~DIAT() { };

	DIAT(const DIAT& source);

	const dArrayT& GetVector_a(){ return a;};
	const dArrayT& GetVector_b(){ return b;};	
	const dArrayT& GetVector_c(){ return c;};

        const dArrayT& GetLatticeParameters();
        const dArray2DT& GetBasis();
	const dArray2DT& GetAxis();

};

#endif
