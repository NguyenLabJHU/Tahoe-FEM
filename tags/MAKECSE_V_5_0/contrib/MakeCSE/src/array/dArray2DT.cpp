/*
 * File: dArray2DT.h
 */

/*
 * created      : PAK (07/16/96)
 * last modified: PAK (05/23/97)
 */

#include <iostream.h>
#include <iomanip.h>

#include "Constants.h"

#include "dArray2DT.h"
#include "iArrayT.h"
#include "dArrayT.h"

/* constructor */
dArray2DT::dArray2DT(void) { }
dArray2DT::dArray2DT(int majordim, int minordim):
	nArray2DT<double>(majordim, minordim) { }
dArray2DT::dArray2DT(int majordim, int minordim, double* p):
	nArray2DT<double>(majordim, minordim, p) { }
dArray2DT::dArray2DT(const dArray2DT& source):
	nArray2DT<double>(source) { }
