/* $Id: SurfaceShapeT.h,v 1.1.1.1 2001-01-29 08:20:31 paklein Exp $ */
/* created: paklein (11/21/1997)                                          */
/* Class to manage CSE integrals, where the dimension of                  */
/* the field variable is 1 greater than the dimension of the parent       */
/* domain. Jump quantities imply jump between any field variable          */
/* across the CSE.                                                        */
/* NOTE: Class operates in 2 modes depending on the dimension             */
/* of coords which are passed in during construction:                     */
/* (1) coords.NumNodes() == fNumFacetNodes: coords used                   */
/* directly as the facet geometry                                         */
/* (2) coords.NumNodes() == fTotalNodes: the facet geometry               */
/* is assumed to be the average of the coordinates                        */
/* on the upper and lower facets.                                         */

#ifndef _SURFACE_SHAPE_T_H_
#define _SURFACE_SHAPE_T_H_

/* base class */
#include "DomainIntegrationT.h"

/* direct members */
#include "iArray2DT.h"
#include "Array2DT.h"
#include "LocalArrayT.h"

class SurfaceShapeT: public DomainIntegrationT
{
public:

	/* constructors */	
	SurfaceShapeT(GeometryT::CodeT geometry_code, int num_ip, int num_nodes,
		int field_dim, const LocalArrayT& coords);
	SurfaceShapeT(const SurfaceShapeT& link, const LocalArrayT& coords);
		// synchronized during integration and shared parent domain,
		// but different coordinates.

	/* accessors */
	int TotalNodes(void) const;
	int NumFacetNodes(void) const;
	int FieldDim(void) const;

	/* set all local parameters */
	virtual void Initialize(void);

/**** for the current integration point ***/

	/* jump in the nodal values */
	const dArrayT& InterpolateJumpU(const LocalArrayT& nodal);

	/* integration point coordinates */
	const dArrayT& IPCoords(void);

	/* extrapolate integration point values to the nodes
	 *    IPvalues[numvals] : values from a single integration point
	 *    nodalvalues[fNumNodes x numvals] : extrapolated values */
	void Extrapolate(const dArrayT& IPvalues, dArray2DT& nodalvalues);

	/* jump gradient tables:
	 *
	 *     fgrad_d = d delta_i/d u_j	[i] = FieldDim
	 *                              	[j] = NumNodes*FieldDim
*
	 *     fgrad_dTgrad_d = d delta_k/d u_i	d delta_k/d u_j
	 *                          	[i],[j] = NumNodes*FieldDim
	 */
	const dMatrixT& Grad_d(void) const;
	const dMatrixT& Grad_dTGrad_d(void) const;

	/* jacobian of the area transformation using the nodes on the 1st facet */
	double Jacobian(void);
	double Jacobian(dMatrixT& Q);
	double Jacobian(dMatrixT& Q, ArrayT<dMatrixT>& dQ);

/*******************************************/

	/* local node numbers on each facet */
	const iArray2DT& NodesOnFacets(void) const;

private:

	/* configure work space arrays */
	void Construct(void);

	/* set jump vector, i.e., assignment of element nodes to
	 * facets: +1 => facet_2, -1 => facet_1. Displacement jump
	 * is: u_2 - u_1. */
	 void SetJumpVector(iArrayT& jump) const;

	/* local node numbers on each facet */
	void SetNodesOnFacets(iArray2DT& facetnodes);

	/* compute average of the facet coordinates */
	void ComputeFacetCoords(void);

private:

	/* dimensions */
	int fTotalNodes;   // twice the facet nodes
	int fNumFacetNodes;// number nodes on each facet
	int fFieldDim;

	/* surface coordinates */
	const LocalArrayT& fCoords;
	LocalArrayT fFacetCoords;

	/* jump shape functions */
	dArray2DT fjumpNa;

	/* shape function tables */
	ArrayT<dMatrixT> fgrad_d;
	ArrayT<dMatrixT> fgrad_dTgrad_d;
	
	/* shape function derivative tables */
	Array2DT<dMatrixT> fgrad_dd;
	
	/* return value */
	dArrayT fInterp;
	
	/* coordinate transformation */
	dMatrixT fJacobian;
	
	/* surface node numbering */
	iArray2DT fFacetNodes;
	dArray2DT fNodalValues; // used for nodal extrapolation
	
	/* work space */
	dArrayT fu_vec;
	dArrayT fx_vec; // shallow
	
	/* work space for 3D jacobian and derivatives */
	dMatrixT fM1;
	dMatrixT fM2;	
};

/* inlines */
inline int SurfaceShapeT::TotalNodes(void) const { return fTotalNodes; }
inline int SurfaceShapeT::NumFacetNodes(void) const { return fNumFacetNodes; }
inline int SurfaceShapeT::FieldDim(void) const { return fFieldDim; }

inline const dArrayT& SurfaceShapeT::IPCoords(void)
{
	/* compute facet coordinates */
	if (fCurrIP == 0 && fCoords.NumberOfNodes() != fNumFacetNodes)
		ComputeFacetCoords();

	fDomain->Interpolate(fFacetCoords, fInterp, fCurrIP);
	return fInterp;
}

inline const dMatrixT& SurfaceShapeT::Grad_d(void) const
{
	return fgrad_d[fCurrIP];
}

inline const dMatrixT& SurfaceShapeT::Grad_dTGrad_d(void) const
{
	return fgrad_dTgrad_d[fCurrIP];
}

/* jacobian of the area transformation at the current integration
* using the nodes on the specified facet */
inline double SurfaceShapeT::Jacobian(void)
{
	/* compute facet coordinates */
	if (fCurrIP == 0 && fCoords.NumberOfNodes() != fNumFacetNodes)
		ComputeFacetCoords();

	/* Jacobian matrix of the surface transformation */
	fDomain->DomainJacobian(fFacetCoords, fCurrIP, fJacobian);	
	return fDomain->SurfaceJacobian(fJacobian);
}

inline double SurfaceShapeT::Jacobian(dMatrixT& Q)
{
	/* compute facet coordinates */
	if (fCurrIP == 0 && fCoords.NumberOfNodes() != fNumFacetNodes)
		ComputeFacetCoords();

	/* Jacobian matrix of the surface transformation */
	fDomain->DomainJacobian(fFacetCoords, fCurrIP, fJacobian);	
	return fDomain->SurfaceJacobian(fJacobian, Q);
}

/* local node numbers on each facet */
inline const iArray2DT& SurfaceShapeT::NodesOnFacets(void) const
{
	return fFacetNodes;
}

#endif /* _SURFACE_SHAPE_T_H_ */
