/* $Id: SSMatSupportT.cpp,v 1.2.2.1 2002-12-10 17:06:03 paklein Exp $ */
#include "SSMatSupportT.h"
#include "ElementsConfig.h"

#ifdef CONTINUUM_ELEMENT
#include "SmallStrainT.h"
#endif

using namespace Tahoe;

/* constructor */
SSMatSupportT::SSMatSupportT(int nsd, int ndof, int nip):
	StructuralMatSupportT(nsd, ndof, nip),
	fStrain_List(NULL),
	fStrain_last_List(NULL),
	fSmallStrain(NULL)
{

}
 
/* destructor */
SSMatSupportT::~SSMatSupportT(void)
{

}

/* set source for the strain */
void SSMatSupportT::SetLinearStrain(const ArrayT<dSymMatrixT>* strain_List)
{
//NOTE: cannot do dimension checks because source is not initialized
//      when this is configured 
#if 0
	/* checks */
	if (!strain_List) throw ExceptionT::kGeneralFail;
	if (strain_List->Length() != NumIP()) throw ExceptionT::kSizeMismatch;
	if (NumIP() > 0 && (*strain_List)[0].Rows() != NumSD()) 
		    throw ExceptionT::kSizeMismatch;
#endif

	/* keep pointer */
	fStrain_List = strain_List;
}

/** set source for the strain from the end of the previous time step */
void SSMatSupportT::SetLinearStrain_last(const ArrayT<dSymMatrixT>* strain_last_List)
{
//NOTE: cannot do dimension checks because source is not initialized
//      when this is configured 
#if 0
	/* checks */
	if (!strain_last_List) throw ExceptionT::kGeneralFail;
	if (strain_last_List->Length() != NumIP()) throw ExceptionT::kSizeMismatch;
	if (NumIP() > 0 && (*strain_last_List)[0].Rows() != NumSD()) 
		    throw ExceptionT::kSizeMismatch;
#endif

	/* keep pointer */
	fStrain_last_List = strain_last_List;
}

/* set the element group pointer */
void SSMatSupportT::SetContinuumElement(const ContinuumElementT* p)
{
	/* inherited */
	StructuralMatSupportT::SetContinuumElement(p);

#ifdef CONTINUUM_ELEMENT
	/* cast to small strain pointer */
	fSmallStrain = dynamic_cast<const SmallStrainT*>(p);
#endif
}
