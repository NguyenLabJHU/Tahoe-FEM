/* $Id: D2FSMatSupportT.cpp,v 1.3 2002-11-15 15:42:04 paklein Exp $ */
#include "D2FSMatSupportT.h"
#include "D2MeshFreeFDElasticT.h"

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

	/* cast to finite strain pointer */
	fD2MeshFreeFDElastic = dynamic_cast<const D2MeshFreeFDElasticT*>(p);
}
