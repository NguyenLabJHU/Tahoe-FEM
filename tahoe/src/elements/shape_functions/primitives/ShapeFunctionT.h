/* $Id: ShapeFunctionT.h,v 1.1.1.1 2001-01-29 08:20:31 paklein Exp $ */
/* created: paklein (06/26/1996)                                          */
/* interface for element shape functions. controls domain representation  */
/* and field representation and spatial derivatives. integration control  */
/* is inherited.                                                          */

#ifndef _SHAPE_FUNCTION_T_H_
#define _SHAPE_FUNCTION_T_H_

/* base class */
#include "DomainIntegrationT.h"

/* direct members */
#include "LocalArrayT.h" // needed for domain

class ShapeFunctionT: public DomainIntegrationT
{
public:

	/* strain-displacement options */
	enum StrainOptionT {kStandardB = 0,
	                  kMeanDilBbar = 1};

/* constructors */
	ShapeFunctionT(GeometryT::CodeT geometry_code, int numIP, const LocalArrayT& coords,
		int B_option);
	ShapeFunctionT(const ShapeFunctionT& link, const LocalArrayT& coords);
		// synchronized during integration and shared parent domain,
		// but different coordinates.

	/* accessors */
//	int NumDOF(void) const;	

	/* type of the domain coordinates */
	LocalArrayT::TypeT DomainCoordType(void) const;

	/* compute global shape derivatives */ 	
	virtual void SetDerivatives(void);

	/* data for all integration points at once */
	const double* IPDets(void) const; // d(fCoords) = j d(parent domain)

/**** for the current integration point ***/
	double IPDet(void) const; // d(fCoords) = j d(parent domain)

	void IPCoords(dArrayT& coordinates) const;
	void InterpolateU(const LocalArrayT& nodal, dArrayT& u) const;
	
	const double* IPShapeX(void) const;	// at domain nodes
	const double* IPShapeU(void) const;	// at field nodes

	void GradU(const LocalArrayT& nodal, dMatrixT& grad_U) const;

	/* extrapolate integration point values to the nodes
	 *    IPvalues[numvals] : values from a single integration point
	 *    nodalvalues[fNumNodes x numvals] : extrapolated values */
	void Extrapolate(const dArrayT& IPvalues, dArray2DT& nodalvalues) const;

	/* convert shape function derivatives by applying a chain rule
	 * transformation:
	 *
	 *      d Na / d x_i = (d Na / d X_J) (d X_J/d x_i)
	 */
	void SetChainRule(const dMatrixT& changeofvar, dArray2DT& derivatives);

	/* strain displacement matrix */
	void B(dMatrixT& B_matrix) const;
	void B(const dArray2DT& derivatives, dMatrixT& B_matrix) const; // 0 =>, Hughes (2.8.20)

	/* shape function gradients matrix */
	void GradNa(dMatrixT& grad_Na) const; // Hughes (4.90)
	void GradNa(const dArray2DT& derivatives, dMatrixT& grad_Na) const; // Hughes (4.90)

	//TEMP - need better way to support different types of
	//       "strain"-"displacement" matrices
	void B_q(dMatrixT& B_matrix) const;
	void B_q(const dArray2DT& derivatives, dMatrixT& B_matrix) const;

/*******************************************/

	/* print the shape function values to the output stream */
	virtual void Print(ostream& out) const;

protected:

	/* set Grad_x matrix - valid until next TopIP/NextIP loop */
	void SetGrad_x(const dArray2DT& Grad_x);

	/* replace field shape function (for non-isoparametric) */
	void SetUShapeFunctions(const dArray2DT& NaU, const ArrayT<dArray2DT>& DNaU);	

	/* access to the (geometry) shape function derivatives */
	const ArrayT<dArray2DT>& DNaX(void) const;

private:

	/* configure work space arrays - initializes shape function to be
	 * isoparametric */
	void Construct(void);
	
	/* hide access to DomainIntegrationT function */
	const double* IPShape(void) const;

	/* compute mean dilatation, Hughes (4.5.23) */
	void SetMeanDilatation(void);

private:

	/* strain-displacement option */
	int fB_option;
	dArray2DT fB_workspace;

	/* local coordinates */
const LocalArrayT& fCoords;

	/* global shape function derivatives */
	dArrayT	fDet;	         // d(fCoords) = j d(parent domain)
	ArrayT<dArray2DT> fDNaX; // geometry: d(phi_X)/d(fCoords)

	/* field shape functions */
	const dArray2DT*         pNaU;  // phi_U
	const ArrayT<dArray2DT>* pDNaU; // d(phi_U)/d(fCoords)

	/* return values */
	const dArray2DT* fGrad_x_temp;
	
	/* work space */
	dArrayT fv1, fv2;
};

/* inlines */
//inline int ShapeFunctionT::NumDOF(void) const { return fNumDOF; }

/* type of the domain coordinates */
inline LocalArrayT::TypeT ShapeFunctionT::DomainCoordType(void) const
{
	return fCoords.Type();
}

/* data for all integration points at once */
inline const double* ShapeFunctionT::IPDets(void) const
{
	return fDet.Pointer();
}

/************************ for the current integration point *********************/
inline double ShapeFunctionT::IPDet(void) const
{
#if __option(extended_errorcheck)
/* range checking */
if (fCurrIP < 0 || fCurrIP >= fNumIP) throw eOutOfRange;
#endif

return *(fDet.Pointer() + fCurrIP);
}

/* data for the current integration point */
inline const double* ShapeFunctionT::IPShapeX(void) const
{
return fDomain->Shape(fCurrIP);
}

inline const double* ShapeFunctionT::IPShapeU(void) const
{
return (*pNaU)(fCurrIP);
}

inline void ShapeFunctionT::IPCoords(dArrayT& coordinates) const
{
	fDomain->Interpolate(fCoords, coordinates, fCurrIP);
}

/* spatial gradients */
inline void ShapeFunctionT::GradU(const LocalArrayT& nodal,
	dMatrixT& grad_U) const
{
	if (fCurrIP != -1)
		fDomain->Jacobian(nodal, (*pDNaU)[fCurrIP], grad_U);
	else
		fDomain->Jacobian(nodal, *fGrad_x_temp, grad_U);
}

inline void ShapeFunctionT::B(dMatrixT& B_matrix) const
{
	B((*pDNaU)[fCurrIP], B_matrix);
}

inline void ShapeFunctionT::GradNa(dMatrixT& grad_Na) const
{
	GradNa((*pDNaU)[fCurrIP], grad_Na);
}

inline void ShapeFunctionT::B_q(dMatrixT& B_matrix) const
{
	B_q((*pDNaU)[fCurrIP], B_matrix);
}

/********************************************************************************/

/* set Grad_x matrix - valid until next TopIP/NextIP loop */
inline void ShapeFunctionT::SetGrad_x(const dArray2DT& Grad_x)
{
	fGrad_x_temp = &Grad_x;

	/* disable IP counter */
	fCurrIP = -1;
}

/* replace field shape function (for non-isoparametric) */
inline void ShapeFunctionT::SetUShapeFunctions(const dArray2DT& NaU,
	const ArrayT<dArray2DT>& DNaU)
{
#if __option(extended_errorcheck)
	if (DNaU.Length() != NumIP()) throw eSizeMismatch;
#endif

	/* set pointers to external data */
	pNaU  = &NaU;
	pDNaU = &DNaU;
}

/* access to the (geometry) shape function derivatives */
inline const ArrayT<dArray2DT>& ShapeFunctionT::DNaX(void) const { return fDNaX; }

#endif /* _SHAPE_FUNCTION_T_H_ */
