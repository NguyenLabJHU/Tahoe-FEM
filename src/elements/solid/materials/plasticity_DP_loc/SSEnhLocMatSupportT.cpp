/* $Id: SSEnhLocMatSupportT.cpp,v 1.1 2005-07-19 18:03:41 raregue Exp $ */

#include "SSEnhLocMatSupportT.h"

#ifdef ENHANCED_STRAIN_LOC_DEV
#include "SmallStrainEnhLocT.h"
#endif

using namespace Tahoe;

/* constructor */
#ifdef ENHANCED_STRAIN_LOC_DEV
SSEnhLocMatSupportT::SSEnhLocMatSupportT(int ndof, int nip):
	SSMatSupportT(ndof, nip),
	fSmallStrainEnhLoc(NULL)
{

}
#else
SSEnhLocMatSupportT::SSEnhLocMatSupportT(int ndof, int nip):
	SSMatSupportT(ndof, nip)
{

}
#endif
 
/* destructor */
SSEnhLocMatSupportT::~SSEnhLocMatSupportT(void)
{

}


/* set source for the stress */
void SSEnhLocMatSupportT::SetElementStress(const ArrayT<dSymMatrixT>* stress_List)
{
	/* keep pointer */
	fStress_List = stress_List;
}

/* set source for the LocFlag */
void SSEnhLocMatSupportT::SetElementLocFlag(const int* loc_flag)
{
	/* keep pointer */
	fLocFlag = loc_flag;
}

/* set the element group pointer */
void SSEnhLocMatSupportT::SetContinuumElement(const ContinuumElementT* p)
{
	/* inherited */
	SSMatSupportT::SetContinuumElement(p);

	/* cast to small strain embedded discontinuity pointer */
	#ifdef ENHANCED_STRAIN_LOC_DEV
	fSmallStrainEnhLoc = TB_DYNAMIC_CAST(const SmallStrainEnhLocT*, p);
	#endif
}
