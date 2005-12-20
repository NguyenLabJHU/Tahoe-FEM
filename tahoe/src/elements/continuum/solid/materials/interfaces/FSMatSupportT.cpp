/* $Id: FSMatSupportT.cpp,v 1.7 2005-12-20 17:26:37 tdnguye Exp $ */
#include "FSMatSupportT.h"
#include "ElementsConfig.h"

#ifdef CONTINUUM_ELEMENT
#include "FiniteStrainT.h"
#endif

using namespace Tahoe;

/* constructor */
FSMatSupportT::FSMatSupportT(int ndof, int nip):
	SolidMatSupportT(ndof, nip),
	fF_List(NULL),
	fF_last_List(NULL),
	fFiniteStrain(NULL)
{

}

/* interpolate the given field to the current integration point */
bool FSMatSupportT::Interpolate_current(const LocalArrayT& u, dArrayT& u_ip) const
{
#ifdef CONTINUUM_ELEMENT
	if (!fFiniteStrain) 
	{
		u_ip = 0.0;
		return false;
	}
	else
	{
        const ShapeFunctionT& CurrShapes = fFiniteStrain->CurrShapeFunction();
        if(&CurrShapes)
            fFiniteStrain->IP_Interpolate_current(u, u_ip);
        else
            fFiniteStrain->IP_Interpolate(u, u_ip);
		return true;
	}
#else
#pragma unused(u)
	u_ip = 0.0;
	return false;
#endif
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
#ifdef CONTINUUM_ELEMENT
	if (fFiniteStrain) {
		fFiniteStrain->ComputeGradient(u, grad_u);
		return true;
	}
	else
		return false;
#else
#pragma unused(u)
#pragma unused(grad_u)
	return false;
#endif
}

/* compute field gradients with respect to current coordinates */
bool FSMatSupportT::ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u, int ip) const
{
#ifdef CONTINUUM_ELEMENT
	if (fFiniteStrain) {
		fFiniteStrain->ComputeGradient(u, grad_u, ip);
		return true;
	}
	else
		return false;
#else
#pragma unused(u)
#pragma unused(grad_u)
#pragma unused(ip)
	return false;
#endif
}

/* compute field gradients with respect to reference coordinates */	
bool FSMatSupportT::ComputeGradient_reference(const LocalArrayT& u, dMatrixT& grad_u) const
{
#ifdef CONTINUUM_ELEMENT
	if (fFiniteStrain) {
		fFiniteStrain->ComputeGradient_reference(u, grad_u);
		return true;
	}
	else
		return false;
#else
#pragma unused(u)
#pragma unused(grad_u)
	return false;
#endif
}

/* compute field gradients with respect to reference coordinates */
bool FSMatSupportT::ComputeGradient_reference(const LocalArrayT& u, dMatrixT& grad_u, int ip) const
{	
#ifdef CONTINUUM_ELEMENT
	if (fFiniteStrain) {
		fFiniteStrain->ComputeGradient_reference(u, grad_u, ip);
		return true;
	}
	else
		return false;
#else
#pragma unused(u)
#pragma unused(grad_u)
#pragma unused(ip)
	return false;
#endif
}

/* set the element group pointer */
void FSMatSupportT::SetContinuumElement(const ContinuumElementT* p)
{
	/* inherited */
	SolidMatSupportT::SetContinuumElement(p);

#ifdef CONTINUUM_ELEMENT
	/* cast to finite strain pointer */
	fFiniteStrain = TB_DYNAMIC_CAST(const FiniteStrainT*, p);
#endif
}
