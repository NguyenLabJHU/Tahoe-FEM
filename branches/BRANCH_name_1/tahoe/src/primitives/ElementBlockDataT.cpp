/* $Id: ElementBlockDataT.cpp,v 1.3.4.1 2002-06-27 18:04:01 cjkimme Exp $ */

#include "ElementBlockDataT.h"

/* copy behavior for arrays FBC_CardT's */

using namespace Tahoe;

const bool ArrayT<ElementBlockDataT*>::fByteCopy = true;
const bool ArrayT<ElementBlockDataT>::fByteCopy = false;

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
