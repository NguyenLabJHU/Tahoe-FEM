/* $Id: SSMatSupportT.cpp,v 1.1.2.2 2002-10-30 09:18:11 paklein Exp $ */
#include "SSMatSupportT.h"
#include "SmallStrainT.h"

using namespace Tahoe;

/* constructor */
SSMatSupportT::SSMatSupportT(int nsd, int ndof, int nip):
	MaterialSupportT(nsd, ndof, nip),
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
	/* checks */
	if (!strain_List) throw ExceptionT::kGeneralFail;
	if (strain_List->Length() != NumIP()) throw ExceptionT::kSizeMismatch;
	if (NumIP() > 0 && (*strain_List)[0].Rows() != NumSD()) 
		    throw ExceptionT::kSizeMismatch;
	
	/* keep pointer */
	fStrain_List = strain_List;
}

/** set source for the strain from the end of the previous time step */
void SSMatSupportT::SetLinearStrain_last(const ArrayT<dSymMatrixT>* strain_last_List)
{
	/* checks */
	if (!strain_last_List) throw ExceptionT::kGeneralFail;
	if (strain_last_List->Length() != NumIP()) throw ExceptionT::kSizeMismatch;
	if (NumIP() > 0 && (*strain_last_List)[0].Rows() != NumSD()) 
		    throw ExceptionT::kSizeMismatch;
	
	/* keep pointer */
	fStrain_last_List = strain_last_List;
}

/* set the element group pointer */
void SSMatSupportT::SetContinuumElement(const ContinuumElementT* p)
{
	/* inherited */
	MaterialSupportT::SetContinuumElement(p);

	/* cast to small strain pointer */
	fSmallStrain = dynamic_cast<const SmallStrainT*>(p);
}

/* return a pointer the specified local array */
const LocalArrayT* SSMatSupportT::LocalArray(LocalArrayT::TypeT t) const
{
	/* quick exit to inherited */
	if (!fSmallStrain) return MaterialSupportT::LocalArray(t);

	switch (t)
	{
		case LocalArrayT::kLastDisp:
			return &(fSmallStrain->LastDisplacements());
	
		case LocalArrayT::kVel:
			return &(fSmallStrain->Velocities());

		case LocalArrayT::kAcc:
			return &(fSmallStrain->Accelerations());

		default:
			/* inherited */
			return MaterialSupportT::LocalArray(t);
	}
}
