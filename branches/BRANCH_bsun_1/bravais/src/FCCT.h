// DEVELOPMENT
#ifndef _FCC_T_H_
#define _FCC_T_H_

#include "ifstreamT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "CrystalLatticeT.h"

using namespace Tahoe;

class FCCT : public CrystalLatticeT 
{
public:
	FCCT(int nlsd,int nuca,dArrayT alat,
	     int which_rot,dArray2DT mat_rot,
	     double angle);

	~FCCT() { };

	FCCT(const FCCT& source);

	const dArrayT& GetVector_a(){ return a;};
	const dArrayT& GetVector_b(){ return b;};	
	const dArrayT& GetVector_c(){ return c;};

        const dArrayT& GetLatticeParameters();
        const dArray2DT& GetBasis();
	const dArray2DT& GetAxis();
};

#endif
