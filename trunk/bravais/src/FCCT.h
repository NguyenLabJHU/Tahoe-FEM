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
	FCCT(int nlsd,int nuca,double alat,dArrayT vec_rot);

	~FCCT() { };

	FCCT(const FCCT& source);

        const dArrayT& GetLatticeParameters();
        const dArray2DT& GetBasis();
};

#endif
