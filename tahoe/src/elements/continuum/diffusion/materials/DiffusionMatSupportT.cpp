/* $Id: DiffusionMatSupportT.cpp,v 1.3 2003-01-27 07:00:28 paklein Exp $ */
#include "DiffusionMatSupportT.h"
#include "ElementsConfig.h"

#ifdef CONTINUUM_ELEMENT
#include "DiffusionT.h"
#endif

using namespace Tahoe;

/* constructor */
DiffusionMatSupportT::DiffusionMatSupportT(int nsd, int ndof, int nip):
	MaterialSupportT(nsd, ndof, nip),
	fGradient_list(NULL),
	fDiffusion(NULL)
{

}

/* set the source for the gradient information */
void DiffusionMatSupportT::SetGradient(const ArrayT<dArrayT>* gradient_list)
{
//NOTE: cannot do dimension checks because source is not initialized
//      when this is configured 
#if 0
	if (!gradient_list ||
	     gradient_list->Length() != NumIP())
	{
		cout << "\n DiffusionMatSupportT::SetGradient: inconsistent gradient source" 
		     << endl;
		throw ExceptionT::kGeneralFail;
	}
#endif

	fGradient_list = gradient_list;
}

/* set the element group pointer */
void DiffusionMatSupportT::SetContinuumElement(const ContinuumElementT* p)
{
	/* inherited */
	MaterialSupportT::SetContinuumElement(p);

#ifdef CONTINUUM_ELEMENT
	/* cast to small strain pointer */
	fDiffusion = dynamic_cast<const DiffusionT*>(p);
#endif
}
