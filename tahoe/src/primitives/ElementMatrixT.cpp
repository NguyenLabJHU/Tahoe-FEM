/* $Id: ElementMatrixT.cpp,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (03/17/1998)                                          */

#include "ElementMatrixT.h"
#include <iostream.h>
#include <iomanip.h>
#include "Constants.h"

/* constructors */
ElementMatrixT::ElementMatrixT(FormatT format):
	fFormat(format)
{

}

ElementMatrixT::ElementMatrixT(int numrows, int numcols, FormatT format):
	dMatrixT(numrows, numcols),
	fFormat(format)
{

}	

ElementMatrixT::ElementMatrixT(int squaredim, FormatT format):
	dMatrixT(squaredim),
	fFormat(format)
{

}

ElementMatrixT::ElementMatrixT(const ElementMatrixT& source):
	dMatrixT(source),
	fFormat(source.fFormat)
{

}

/* take symmetric with values in upper only and copy to full */
/* NOTE: error to call for kNonSymmetric matrices            */
/* NOTE: not strictly const, but in essence                  */
void ElementMatrixT::CopySymmetric(void) const
{
#if __option (extended_errorcheck)
	if (fFormat == kNonSymmetric) throw eGeneralFail;
#endif

	if (fFormat == kDiagonal)
		return;
	else
	{
		/* cast away const-ness */
		dMatrixT* localthis = (dMatrixT*) this;
		localthis->CopySymmetric();
	}
}
