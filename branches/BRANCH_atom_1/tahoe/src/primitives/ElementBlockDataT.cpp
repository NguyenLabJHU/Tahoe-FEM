/* $Id: ElementBlockDataT.cpp,v 1.5 2002-07-05 17:18:39 paklein Exp $ */

#include "ElementBlockDataT.h"

/* copy behavior for arrays FBC_CardT's */

using namespace Tahoe;

namespace Tahoe {
const bool ArrayT<ElementBlockDataT*>::fByteCopy = true;
const bool ArrayT<ElementBlockDataT>::fByteCopy = false;
} /* namespace Tahoe */

/* constructor */
ElementBlockDataT::ElementBlockDataT(void):
	fStartNum(-1),
	fDimension(-1),
	fMaterial(-1)
{

}

/* copy constructor */
ElementBlockDataT::ElementBlockDataT(ElementBlockDataT& source)
{
	operator=(source);
}

/* set field */
void ElementBlockDataT::Set(const StringT& ID, int start, int dim, int material)
{
	fID = ID;
	fStartNum = start;
	fDimension = dim;
	fMaterial = material;
}
