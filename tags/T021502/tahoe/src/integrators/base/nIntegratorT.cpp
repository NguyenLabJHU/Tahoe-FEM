/* $Id: nIntegratorT.cpp,v 1.2 2001-08-27 17:12:12 paklein Exp $ */
/* created: paklein (10/14/1996) */

#include "nIntegratorT.h"
#include <iostream.h>
#include "ExceptionCodes.h"

/* constructor */
nIntegratorT::nIntegratorT(int order)
{
	if (order < 0) throw eGeneralFail;
	fU.Allocate(order + 1);
	fU = NULL;
}

/* destructor */
nIntegratorT::~nIntegratorT(void) { }

/* register field arrays */
void nIntegratorT::RegisterField(dArray2DT& field, int order)
{
	/* check */
	if (order < 0 || order > fU.Length())
	{
		cout << "\n nIntegratorT::RegisterField: derivative is out of range {0,"
		     << fU.Length() - 1 << "}:" << order << endl;
		throw eOutOfRange;
	}
	
	/* store pointer */
	fU[order] = &field;
}
