// DEVELOPMENT
/* $Id: PeriodicTableT.h,v 1.4 2002-11-14 01:47:33 saubry Exp $ */

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

	ArrayT <PerTabEntryT> PT;

	PeriodicTableT();
	~PeriodicTableT(){};

	void Initialize();
	PerTabEntryT operator[] (const char * s);
	PerTabEntryT operator[] (const int m);
};

#endif
