/* $Id: DiffusionMatSupportT.cpp,v 1.6 2003-12-28 08:23:16 paklein Exp $ */
#include "DiffusionMatSupportT.h"
#include "ElementsConfig.h"

#ifdef CONTINUUM_ELEMENT
#include "DiffusionElementT.h"
#endif

using namespace Tahoe;

/* constructor */
DiffusionMatSupportT::DiffusionMatSupportT(int nsd, int ndof, int nip):
	MaterialSupportT(nsd, ndof, nip),
	fField_list(NULL),
	fGradient_list(NULL),
	fDiffusion(NULL)
{

}

/* set the source for the gradient information */
void DiffusionMatSupportT::SetField(const dArrayT* field_list)
{
	fField_list = field_list;
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
	fDiffusion = TB_DYNAMIC_CAST(const DiffusionElementT*, p);
#endif
}
