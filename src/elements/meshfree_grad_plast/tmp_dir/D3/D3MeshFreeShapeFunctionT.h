/* $Id: D3MeshFreeShapeFunctionT.h,v 1.1 2004-09-28 23:13:59 raregue Exp $ */
/* created: paklein (10/23/1999) */
#ifndef _D3_MF_SHAPE_T_H_
#define _D3_MF_SHAPE_T_H_

/* base class */
#include "D2MeshFreeShapeFunctionT.h"

namespace Tahoe {

/* forward declarations */
class D3MeshFreeSupportT;

class D3MeshFreeShapeFunctionT: public D2MeshFreeShapeFunctionT
{
public:

/* constructors */
	D3MeshFreeShapeFunctionT(GeometryT::CodeT geometry_code, int numIP,
		const LocalArrayT& coords, const dArray2DT& all_coords,
		const iArray2DT& connects, const iArrayT& nongridnodes,
		const int& currelement, const ParameterListT& mf_support_params);

	/* compute global shape derivatives */ 	
	virtual void SetDerivatives(void);
	int SetDerivativesAt(const dArrayT& x, AutoArrayT<int>& nodes);
		// returns 0 if MLS fit fails

	/* 3rd order shape function gradients matrix */
	void D3GradNa(dMatrixT& D3_grad_Na) const;

	/* 3rd order spatial gradients */
	void GradGradGradU(const LocalArrayT& nodal, dMatrixT& gradgradgrad_U) const;

	/* 3rd derivatives of shape functions at IP */
	const dArray2DT& DDDerivatives_U(int ip) const { return fDDDNaU[ip]; };
	const dArray2DT& DDDerivatives_U(void) const { return fDDDNaU[fCurrIP]; };

	/* reconstruct displacement field and all derivatives */
	void NodalField(const dArray2DT& DOF, dArray2DT& field, dArray2DT& Dfield,
		dArray2DT& DDfield, dArray2DT& DDDfield, iArrayT& nodes);

protected:

	/* meshfree database support */
	D3MeshFreeSupportT* fD3MFSupport;

	ArrayT<dArray2DT> fDDDNaU;
	
	/* work space for blended shape functions */
	ArrayT<dArray2DT> fDDDNa_tmp;
};

/* inlines */

/* spatial gradients */
inline void D3MeshFreeShapeFunctionT::GradGradGradU(const LocalArrayT& nodal,
	dMatrixT& gradgradgrad_U) const
{
	//fDomain->??(nodal, fDDDNaU[fCurrIP], gradgradgrad_U);	
}

} // namespace Tahoe 
#endif /* _D3_MF_SHAPE_T_H_ */
