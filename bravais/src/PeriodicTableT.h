/* $Id: PeriodicTableT.h,v 1.1 2002-03-13 02:25:55 jzimmer Exp $ */

#ifndef _PERIODIC_TABLE_T_H_
#define _PERIODIC_TABLE_T_H_

#include "dArrayT.h"
#include "dArray2DT.h"
#include "StringT.h"
#include "ArrayT.h"
#include "Array2DT.h"
#include "PerTabEntryT.h"

class PeriodicTableT {
public:
	ArrayT <PerTabEntryT> PT(110);	
	
	PeriodicTableT();
	~PeriodicTableT();

	void Initialize();
	PerTabEntryT operator[] (const char * s);
};

#endif
