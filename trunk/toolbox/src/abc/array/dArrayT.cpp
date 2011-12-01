/* $Id: dArrayT.cpp,v 1.12 2011-12-01 20:25:15 bcyansfn Exp $ */
/* created: paklein (08/11/1996) */
#include "dArrayT.h"
#include <iostream>
#include <iomanip>
#include <cmath>

#include "toolboxConstants.h"

using namespace Tahoe;

namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<dArrayT*>::fByteCopy = true; 
DEFINE_TEMPLATE_STATIC const bool ArrayT<dArrayT>::fByteCopy = false; 
} /* namespace Tahoe */

/* constructor */
dArrayT::dArrayT(void) { }
dArrayT::dArrayT(int length): nArrayT<double>(length) { }
dArrayT::dArrayT(int length, const double* p): nArrayT<double>(length,p) { }
dArrayT::dArrayT(const dArrayT& source): nArrayT<double>(source) { }

/* L2 norm of the vector */
double dArrayT::Magnitude(void) const
{
	int length = Length();
	const double* p = Pointer();
	if (length == 2)
		return sqrt(p[0]*p[0] + p[1]*p[1]);
	else if (length == 3)
		return sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
	else if (length == 1)
		return fabs(p[0]);
	else
	{
		register double magsqr = 0.0;
		for (int i = 0; i < Length(); i++)
		{
			magsqr += (*p)*(*p);
			p++;
		}
		return sqrt(magsqr);
	}
}
