/* $Id: ElementBlockDataT.cpp,v 1.6.10.1 2003-10-27 22:22:52 paklein Exp $ */
#include "ElementBlockDataT.h"

using namespace Tahoe;

/* copy behavior for arrays FBC_CardT's */
namespace Tahoe {
template<> const bool ArrayT<ElementBlockDataT*>::fByteCopy = true;
template<> const bool ArrayT<ElementBlockDataT>::fByteCopy = false;
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
