/* $Id: IC_CardT.cpp,v 1.14.2.4 2004-07-13 16:48:23 paklein Exp $ */
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
IC_CardT::IC_CardT(void): fnode(-1), fdof(-1), fvalue(0.0)			
{
	//initialize to inappropriate values
}

void IC_CardT::SetValues(int node, int dof, int order, double value)
{
	/* set */
	fnode  = node;
	fdof   = dof;
	forder = order;
	fvalue = value;			

	/* check */
	if (order < 0) throw ExceptionT::kBadInputValue;
}
