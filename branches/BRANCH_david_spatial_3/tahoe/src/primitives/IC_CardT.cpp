/* $Id: IC_CardT.cpp,v 1.15.30.1 2005-07-25 02:37:24 paklein Exp $ */
/* created: paklein (07/16/1997) */
#include "IC_CardT.h"
#include "ArrayT.h"

using namespace Tahoe;

/* copy behavior for arrays IC_CardT's */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<IC_CardT*>::fByteCopy = true;
DEFINE_TEMPLATE_STATIC const bool ArrayT<IC_CardT>::fByteCopy = false;
} /* namespace Tahoe */

/* default constructor */
IC_CardT::IC_CardT(void): 
	fnode(-1), 
	fdof(-1), 
	fvalue(0.0),
	fmode(kUndefined)
{
}

void IC_CardT::SetValues(int node, int dof, int order, double value)
{
	/* set */
	fnode  = node;
	fdof   = dof;
	forder = order;
	fvalue = value;
	fmode = kNode;
	fID.Clear();

	/* check */
	if (order < 0) ExceptionT::BadInputValue("IC_CardT::SetValues");
}

void IC_CardT::SetValues(const StringT& ID, int dof, int order, double value)
{
	/* set */
	fnode  = -1;
	fdof   = dof;
	forder = order;
	fvalue = value;
	fmode  = kSet;
	fID    = ID;

	/* check */
	if (order < 0) ExceptionT::BadInputValue("IC_CardT::SetValues");
}
