#ifndef _DIA_T_H_
#define _DIA_T_H_

#include "StringT.h"
#include "ExceptionCodes.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "CrystalLatticeT.h"
#include "ifstreamT.h"

using namespace Tahoe;

class DIAT : public CrystalLatticeT 
{
public:
	DIAT(int nlsd,int nuca,double alat);

	~DIAT() { };

	DIAT(const DIAT& source);

        const dArrayT& GetLatticeParameters();
        const dArray2DT& GetBasis();
};

#endif
