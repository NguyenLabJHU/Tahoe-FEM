/* $Id: SimoShapeFunctionT.h,v 1.4 2001-07-20 00:58:39 paklein Exp $ */

#ifndef _SIMO_SHAPE_FUNCTION_T_H_
#define _SIMO_SHAPE_FUNCTION_T_H_

/* base class */
#include "ShapeFunctionT.h"

/** enhanced modes shape functions. for use in conjunction with element
 * formulation due to Simo, Armero, and Taylor, CMAME \b 110, 359-386, 1993 
 * \note equation numbers refer to the paper. */
class SimoShapeFunctionT: public ShapeFunctionT
{
public:

	/** constructor */
	SimoShapeFunctionT(GeometryT::CodeT geometry_code, int numIP, 
		const LocalArrayT& coords, const LocalArrayT& element_modes);

	/** compute global shape derivatives */ 	
	virtual void SetDerivatives(void);

	/** strain displacement matrix associated with the enhanced modes (4.11)
	 * at the current integration point */
	void B_enhanced(dMatrixT& B_matrix) const;

	/** shape function gradients matrix associated with the enhanced 
	 * modes, as in (4.18) at the current integration point */
	void GradNa_enhanced(dMatrixT& grad_Na) const;

	/** enhanced part of the field gradients at the current integration point.
	 * This part of the deformation gradient is associated only with the
	 * bubble enhancement.
	 * \param nodal array of nodal values: [nnd] x [nu]
	 * \param grad_U field gradient matrix: [nu] x [nsd] */
	void GradU_enhanced(const LocalArrayT& nodal, dMatrixT& grad_U) const;

	/** enhanced part of the field gradients at the specified integration point. 
	 * This part of the deformation gradient is associated only with the
	 * bubble enhancement.
	 * \param nodal array of nodal values: [nnd] x [nu]
	 * \param grad_U field gradient matrix: [nu] x [nsd] 
	 * \param IPnumber integration point number */
	void GradU_enhanced(const LocalArrayT& nodal, dMatrixT& grad_U, int IPnumber) const;

	/** write shape function values to the output stream */
	virtual void Print(ostream& out) const;

private:

	/** element modes */
	const LocalArrayT& fElementModes;
	bool fHas3DIncompressibleMode;

	/** gradients of the enhanced bubble modes */
	ArrayT<dArray2DT> fDNaX_bubble;

	/** gradients of the enhanced incompressible modes (3D only) */
	dArray2DT fDNaX_inc;
};

/* inlines */

/* strain displacement matrix associated with the enhanced modes (4.11) */
inline void SimoShapeFunctionT::B_enhanced(dMatrixT& B_matrix) const
{
	/* inherited */
	B(fDNaX_bubble[fCurrIP], B_matrix);
}

/* shape function gradients matrix */
inline void SimoShapeFunctionT::GradNa_enhanced(dMatrixT& grad_Na) const
{
	/* inherited */
	GradNa(fDNaX_bubble[fCurrIP], grad_Na);
}

/* enhanced part of the field gradients at the current integration point */
inline void SimoShapeFunctionT::GradU_enhanced(const LocalArrayT& nodal, 
	dMatrixT& grad_U) const
{
	/* computed  by parent domain */
	fDomain->Jacobian(nodal, fDNaX_bubble[fCurrIP], grad_U);
}

inline void SimoShapeFunctionT::GradU_enhanced(const LocalArrayT& nodal, 
	dMatrixT& grad_U, int IPnumber) const
{
	/* computed  by parent domain */
	fDomain->Jacobian(nodal, fDNaX_bubble[IPnumber], grad_U);
}

#endif /* _SIMO_SHAPE_FUNCTION_T_H_ */
