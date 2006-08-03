/* $Id: FSFiberMatT.cpp,v 1.1 2006-08-03 01:10:41 thao Exp $ */
/* created: paklein (06/09/1997) */
#include "FSFiberMatT.h"
#include "FSFiberMatSupportT.h"
#include "iArray2DT.h"

using namespace Tahoe;

/* constructor */
FSFiberMatT::FSFiberMatT(void):
	ParameterInterfaceT("fiber_composite_material"),
	fFSFiberMatSupport(NULL)
{

}

/* set the material support or pass NULL to clear */
void FSFiberMatT::SetFSFiberMatSupport(const FSFiberMatSupportT* support)
{
	/* set inherited material support */
	FSSolidMatT::SetFSMatSupport(support);

	fFSFiberMatSupport = support;

}