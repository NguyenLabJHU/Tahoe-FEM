/* $Id: D2FSMatSupportT.cpp,v 1.4 2003-01-27 07:00:28 paklein Exp $ */
#include "D2FSMatSupportT.h"
#include "ElementsConfig.h"

#ifdef CONTINUUM_ELEMENT
#include "D2MeshFreeFDElasticT.h"
#endif

using namespace Tahoe;

/* constructor */
D2FSMatSupportT::D2FSMatSupportT(int nsd, int ndof, int nip):
	FSMatSupportT(nsd, ndof, nip),
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
	fD2MeshFreeFDElastic = dynamic_cast<const D2MeshFreeFDElasticT*>(p);
#endif
}
