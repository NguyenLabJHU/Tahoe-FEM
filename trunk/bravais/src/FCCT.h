/* $Id: FCCT.h,v 1.2 2002-07-24 01:15:15 saubry Exp $ */

#ifndef _FCC_T_H_
#define _FCC_T_H_

#include "StringT.h"
#include "ExceptionCodes.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "CrystalLatticeT.h"
#include "ifstreamT.h"

using namespace Tahoe;

class FCCT : public CrystalLatticeT 
{
public:
	FCCT(int nlsd,int nuca,double alat);

	~FCCT() { };

	FCCT(const FCCT& source);

        const dArrayT& GetLatticeParameters();
        const dArray2DT& GetBasis();
};

#endif
