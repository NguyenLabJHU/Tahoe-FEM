/* $Id: nIntegratorT.cpp,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (10/14/1996)                                          */

#include "nIntegratorT.h"
#include <iostream.h>
#include "ExceptionCodes.h"

/* constructor */
nIntegratorT::nIntegratorT(void): fU(NULL) { }

/* destructor */
nIntegratorT::~nIntegratorT(void) { }

/* register field arrays */
void nIntegratorT::SetField(dArray2DT& field, int order)
{
	if (order == 0)
		fU = &field;
	else
	{
		cout << "\n nIntegratorT::SetField: unexpected field order: "
		     << order << endl;
		throw eOutOfRange;
	}
}
