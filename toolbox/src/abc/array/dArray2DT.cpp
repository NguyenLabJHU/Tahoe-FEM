/* $Id: dArray2DT.cpp,v 1.4.2.1 2002-06-27 18:00:45 cjkimme Exp $ */
/* created: paklein (07/16/1996) */

#include "dArray2DT.h"

#include <iostream.h>
#include <iomanip.h>

#include "Constants.h"
#include "iArrayT.h"
#include "dArrayT.h"

using namespace Tahoe;

/* array behavior */
const bool ArrayT<dArray2DT*>::fByteCopy = true;
const bool ArrayT<const dArray2DT*>::fByteCopy = true;
const bool ArrayT<dArray2DT>::fByteCopy  = false;

/* constructor */
dArray2DT::dArray2DT(void) { }
dArray2DT::dArray2DT(int majordim, int minordim):
	nArray2DT<double>(majordim, minordim) { }
dArray2DT::dArray2DT(int majordim, int minordim, double* p):
	nArray2DT<double>(majordim, minordim, p) { }
dArray2DT::dArray2DT(const dArray2DT& source):
	nArray2DT<double>(source) { }
