/* $Id: FSMatSupportT.cpp,v 1.2 2002-11-14 17:06:21 paklein Exp $ */
#include "FSMatSupportT.h"
#include "FiniteStrainT.h"

using namespace Tahoe;

/* constructor */
FSMatSupportT::FSMatSupportT(int nsd, int ndof, int nip):
	StructuralMatSupportT(nsd, ndof, nip),
	fF_List(NULL),
	fF_last_List(NULL),
	fFiniteStrain(NULL)
{

}

/* set source for the deformation gradient */
void FSMatSupportT::SetDeformationGradient(const ArrayT<dMatrixT>* F_List)
{
//NOTE: cannot do dimension checks because source is not initialized
//      when this is configured 
#if 0
	/* checks */
	if (!F_List) throw ExceptionT::kGeneralFail;
	if (F_List->Length() != NumIP()) throw ExceptionT::kSizeMismatch;
	if (NumIP() > 0)
		if ((*F_List)[0].Rows() != NumSD() || /* only check the first one */
		    (*F_List)[0].Cols() != NumSD()) 
		    throw ExceptionT::kSizeMismatch;
#endif
	
	/* keep pointer */
	fF_List = F_List;
}

/* set source for the deformation gradient */
void FSMatSupportT::SetDeformationGradient_last(const ArrayT<dMatrixT>* F_last_List)
{
//NOTE: cannot do dimension checks because source is not initialized
//      when this is configured 
#if 0
	/* checks */
	if (!F_last_List) throw ExceptionT::kGeneralFail;
	if (F_last_List->Length() != NumIP()) throw ExceptionT::kSizeMismatch;
	if (NumIP() > 0)
		if ((*F_last_List)[0].Rows() != NumSD() || /* only check the first one */
		    (*F_last_List)[0].Cols() != NumSD()) 
		    throw ExceptionT::kSizeMismatch;
#endif

	/* keep pointer */
	fF_last_List = F_last_List;
}

/* compute field gradients with respect to current coordinates */
bool FSMatSupportT::ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u) const
{
	if (fFiniteStrain) {
		fFiniteStrain->ComputeGradient(u, grad_u);
		return true;
	}
	else
		return false;
}

/* compute field gradients with respect to current coordinates */
bool FSMatSupportT::ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u, int ip) const
{
	if (fFiniteStrain) {
		fFiniteStrain->ComputeGradient(u, grad_u, ip);
		return true;
	}
	else
		return false;
}

/* compute field gradients with respect to reference coordinates */	
bool FSMatSupportT::ComputeGradient_reference(const LocalArrayT& u, dMatrixT& grad_u) const
{
	if (fFiniteStrain) {
		fFiniteStrain->ComputeGradient_reference(u, grad_u);
		return true;
	}
	else
		return false;
}

/* compute field gradients with respect to reference coordinates */
bool FSMatSupportT::ComputeGradient_reference(const LocalArrayT& u, dMatrixT& grad_u, int ip) const
{	
	if (fFiniteStrain) {
		fFiniteStrain->ComputeGradient_reference(u, grad_u, ip);
		return true;
	}
	else
		return false;
}

/* set the element group pointer */
void FSMatSupportT::SetContinuumElement(const ContinuumElementT* p)
{
	/* inherited */
	StructuralMatSupportT::SetContinuumElement(p);

	/* cast to finite strain pointer */
	fFiniteStrain = dynamic_cast<const FiniteStrainT*>(p);
}
