/* $Id: FSMatSupportT.cpp,v 1.1.2.2 2002-10-30 09:18:11 paklein Exp $ */
#include "FSMatSupportT.h"
#include "FiniteStrainT.h"

using namespace Tahoe;

/* constructor */
FSMatSupportT::FSMatSupportT(int nsd, int ndof, int nip):
	MaterialSupportT(nsd, ndof, nip),
	fF_List(NULL),
	fF_last_List(NULL),
	fFiniteStrain(NULL)
{

}

/* set source for the deformation gradient */
void FSMatSupportT::SetDeformationGradient(const ArrayT<dMatrixT>* F_List)
{
	/* checks */
	if (!F_List) throw ExceptionT::kGeneralFail;
	if (F_List->Length() != NumIP()) throw ExceptionT::kSizeMismatch;
	if (NumIP() > 0)
		if ((*F_List)[0].Rows() != NumSD() || /* only check the first one */
		    (*F_List)[0].Cols() != NumSD()) 
		    throw ExceptionT::kSizeMismatch;
	
	/* keep pointer */
	F_List = fF_List;
}

/* set source for the deformation gradient */
void FSMatSupportT::SetDeformationGradient_last(const ArrayT<dMatrixT>* F_last_List)
{
	/* checks */
	if (!F_last_List) throw ExceptionT::kGeneralFail;
	if (F_last_List->Length() != NumIP()) throw ExceptionT::kSizeMismatch;
	if (NumIP() > 0)
		if ((*F_last_List)[0].Rows() != NumSD() || /* only check the first one */
		    (*F_last_List)[0].Cols() != NumSD()) 
		    throw ExceptionT::kSizeMismatch;
	
	/* keep pointer */
	fF_last_List = F_last_List;
}

/* set the element group pointer */
void FSMatSupportT::SetContinuumElement(const ContinuumElementT* p)
{
	/* inherited */
	MaterialSupportT::SetContinuumElement(p);

	/* cast to finite strain pointer */
	fFiniteStrain = dynamic_cast<const FiniteStrainT*>(p);
}

/* return a pointer the specified local array */
const LocalArrayT* FSMatSupportT::LocalArray(LocalArrayT::TypeT t) const
{
	/* quick exit to inherited */
	if (!fFiniteStrain) return MaterialSupportT::LocalArray(t);

	switch (t)
	{
		case LocalArrayT::kLastDisp:
			return &(fFiniteStrain->LastDisplacements());
	
		case LocalArrayT::kVel:
			return &(fFiniteStrain->Velocities());

		case LocalArrayT::kAcc:
			return &(fFiniteStrain->Accelerations());

		default:
			/* inherited */
			return MaterialSupportT::LocalArray(t);
	}
}

/* nodal temperatures */
const LocalArrayT* FSMatSupportT::Temperatures(void) const
{
	if (!fFiniteStrain)
		return NULL;
	else
		return fFiniteStrain->Temperatures();
}

/* nodal temperatures from the last time step */
const LocalArrayT* FSMatSupportT::LastTemperatures(void) const
{
	if (!fFiniteStrain)
		return NULL;
	else
		return fFiniteStrain->LastTemperatures();
}
