/* $Id: IC_CardT.cpp,v 1.14 2004-06-17 07:14:05 paklein Exp $ */
/* created: paklein (07/16/1997) */
#include "IC_CardT.h"

#include <iostream.h>
#include <iomanip.h>

#include "toolboxConstants.h"
#include "ExceptionT.h"
#include "ifstreamT.h"

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

/* modifiers */
void IC_CardT::SetValues(ifstreamT& in)
{
	/* parameters */
	int    node;
	int    dof;
	int    order;
	double value;

	/* read */
	in >> node >> dof >> order >> value;
	
	/* correct offset */
	if (node > 0) node--;
	dof--;
	
	/* set and echo */
	SetValues(node, dof, order, value);
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

/* I/O */
void IC_CardT::WriteHeader(ostream& out)
{
	int d_width = out.precision() + kDoubleExtra;
	out << setw(kIntWidth)    << "node";
	out << setw(kIntWidth)    << "dof";
	out << setw(kIntWidth)    << "order";
	out << setw(d_width)      << "value" << '\n';
}

void IC_CardT::WriteValues(ostream& out) const
{
	/* output format */
	int d_width = out.precision() + kDoubleExtra;

	/* echo */
	if (fnode == -1)
		out << setw(kIntWidth)    << "(ALL) 0";
	else
		out << setw(kIntWidth)    << fnode + 1;
	out << setw(kIntWidth) << fdof + 1;
	out << setw(kIntWidth) << forder;
	out << setw(d_width)   << fvalue << '\n';
}
