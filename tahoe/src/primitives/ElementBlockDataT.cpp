/* $Id: ElementBlockDataT.cpp,v 1.6 2003-08-14 05:32:30 paklein Exp $ */
#include "ElementBlockDataT.h"

using namespace Tahoe;

/* copy behavior for arrays FBC_CardT's */
namespace Tahoe {
const bool ArrayT<ElementBlockDataT*>::fByteCopy = true;
const bool ArrayT<ElementBlockDataT>::fByteCopy = false;
} /* namespace Tahoe */

/* constructor */
ElementBlockDataT::ElementBlockDataT(void):
	ParameterInterfaceT("element_block"),
	fStartNum(-1),
	fDimension(-1),
	fMaterial(-1)
{

}

/* copy constructor */
ElementBlockDataT::ElementBlockDataT(ElementBlockDataT& source):
	ParameterInterfaceT("element_block")
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

/* describe the parameters needed by the interface */
void ElementBlockDataT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	list.AddParameter(fID, "block_ID");
	list.AddParameter(fMaterial, "material");
}
