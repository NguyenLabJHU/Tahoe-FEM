/* $Id: D2MeshFreeShapeFunctionT.h,v 1.12 2005-12-08 15:32:15 kyonten Exp $ */
/* created: paklein (10/23/1999) */
#ifndef _D2_MF_SHAPE_T_H_
#define _D2_MF_SHAPE_T_H_

/* base class */
#include "MeshFreeShapeFunctionT.h"

namespace Tahoe {

/* forward declarations */
class D2MeshFreeSupportT;

class D2MeshFreeShapeFunctionT: public MeshFreeShapeFunctionT
{
public:

	/* constructors */
	D2MeshFreeShapeFunctionT(GeometryT::CodeT geometry_code, int numIP,
		const LocalArrayT& coords, const dArray2DT& all_coords,
		const iArray2DT& connects, const iArrayT& nongridnodes,
		const int& currelement, const ParameterListT& mf_support_params);

	/** class-dependent initializations */
	virtual void Initialize(void);

	/* compute global shape derivatives */ 	
	virtual void SetDerivatives(void);
	int SetDerivativesAt(const dArrayT& x, AutoArrayT<int>& nodes);
		// returns 0 if MLS fit fails

	/* 2nd order shape function gradients matrix */
	void D2GradNa(dMatrixT& D2_grad_Na) const;

	/* 2nd order spatial gradients */
	void GradGradU(const LocalArrayT& nodal, dMatrixT& gradgrad_U) const;
	void GradGradU(const LocalArrayT& nodal, dMatrixT& gradgrad_U, int ip) const;
	
	/* 2nd derivatives of shape functions at IP */
	const dArray2DT& DDerivatives_U(int ip) const { return fDDNaU[ip]; };
	const dArray2DT& DDerivatives_U(void) const { return fDDNaU[fCurrIP]; };

	/* reconstruct displacement field and all derivatives */
	void NodalField(const dArray2DT& DOF, dArray2DT& field, dArray2DT& Dfield,
		dArray2DT& DDfield, iArrayT& nodes);
		
	/* reference to the support */
	D2MeshFreeSupportT& D2MeshFreeSupport(void) const;

protected:

	/* meshfree database support */
	D2MeshFreeSupportT* fD2MFSupport;

	ArrayT<dArray2DT> fDDNaU;
	
	/* work space for blended shape functions */
	ArrayT<dArray2DT> fDDNa_tmp;
};

/* inlines */
inline D2MeshFreeSupportT& D2MeshFreeShapeFunctionT::D2MeshFreeSupport(void) const
{
	if (!fD2MFSupport) throw ExceptionT::kGeneralFail;
	return *fD2MFSupport;
}

/* spatial gradients */
inline void D2MeshFreeShapeFunctionT::GradGradU(const LocalArrayT& nodal,
	dMatrixT& gradgrad_U) const
{
#if __option(extended_errorcheck)
	if (gradgrad_U.Rows() != nodal.MinorDim() ||
		gradgrad_U.Cols() != fDDNaU[fCurrIP].MajorDim())
    	throw ExceptionT::kSizeMismatch;
#endif
	
	fDomain->Jacobian(nodal, fDDNaU[fCurrIP], gradgrad_U);	
}

inline void D2MeshFreeShapeFunctionT::GradGradU(const LocalArrayT& nodal,
	dMatrixT& gradgrad_U, int ip) const
{
#if __option(extended_errorcheck)
	if (gradgrad_U.Rows() != nodal.MinorDim() ||
		gradgrad_U.Cols() != fDDNaU[ip].MajorDim())
    	throw ExceptionT::kSizeMismatch;
#endif
	
	fDomain->Jacobian(nodal, fDDNaU[ip], gradgrad_U);	
}

} // namespace Tahoe 
#endif /* _D2_MF_SHAPE_T_H_ */
