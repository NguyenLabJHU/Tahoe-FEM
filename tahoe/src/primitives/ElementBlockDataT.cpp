/* $Id: ElementBlockDataT.cpp,v 1.2 2002-02-27 01:32:55 paklein Exp $ */

#include "ElementBlockDataT.h"

/* copy behavior for arrays FBC_CardT's */
template<> const bool ArrayT<ElementBlockDataT*>::fByteCopy = true;
template<> const bool ArrayT<ElementBlockDataT>::fByteCopy = false;

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
