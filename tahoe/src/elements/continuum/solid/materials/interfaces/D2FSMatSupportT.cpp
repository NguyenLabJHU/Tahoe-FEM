/* $Id: D2FSMatSupportT.cpp,v 1.6.18.1 2004-06-14 04:56:33 paklein Exp $ */
#include "D2FSMatSupportT.h"
#include "ElementsConfig.h"

#ifdef CONTINUUM_ELEMENT
#include "D2MeshFreeFSSolidT.h"
#endif

using namespace Tahoe;

/* constructor */
D2FSMatSupportT::D2FSMatSupportT(int ndof, int nip):
	FSMatSupportT(ndof, nip),
	fD2MeshFreeFDElastic(NULL)
{

}

/* set the element group pointer */
void D2FSMatSupportT::SetContinuumElement(const ContinuumElementT* p)
{
	/* inherited */
	FSMatSupportT::SetContinuumElement(p);

#ifdef CONTINUUM_ELEMENT
	/* cast to finite strain pointer */
	fD2MeshFreeFDElastic = TB_DYNAMIC_CAST(const D2MeshFreeFSSolidT*, p);
#endif
}
