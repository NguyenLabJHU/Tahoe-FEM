/* $Id: FCCT.h,v 1.1 2002-03-06 01:55:43 jzimmer Exp $ */

#ifndef _FCC_T_H_
#define _FCC_T_H_

#include "StringT.h"
#include "ExceptionCodes.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "CrystalLatticeT.h"

class ifstreamT;

class FCCT : public CrystalLatticeT {
public:
	FCCT();
	~FCCT() { }

	void SetBasis();
	void SetLatticeParameters(ifstreamT& in);
        const dArrayT& GetLatticeParameters();
        const dArray2DT& GetBasis();
};





#endif
