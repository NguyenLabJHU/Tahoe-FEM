/* $Id: SurfaceShapeT.h,v 1.10 2005-07-29 08:13:02 paklein Exp $ */
/* created: paklein (11/21/1997) */

#ifndef _SURFACE_SHAPE_T_H_
#define _SURFACE_SHAPE_T_H_

/* base class */
#include "DomainIntegrationT.h"

/* direct members */
#include "iArray2DT.h"
#include "Array2DT.h"
#include "LocalArrayT.h"

namespace Tahoe {

/** Class to manage CSE integrals, where the dimension of
 * the field variable is 1 greater than the dimension of the parent
 * domain. Jump quantities imply jump between any field variable
 * across the CSE.
 * \note Class operates in 2 modes depending on the dimension
 * of coords which are passed in during construction:
 * (1) coords.NumNodes() == fNumFacetNodes: coords used
 * directly as the facet geometry
 * (2) coords.NumNodes() == fTotalNodes: the facet geometry
 * is assumed to be the average of the coordinates
 * on the upper and lower facets. */
class SurfaceShapeT: public DomainIntegrationT
{
public:

	/** constructor. 
	 * \param geometry_code geometry of the parent domain
	 * \param num_ip number of integration points 
	 * \param num_nodes total number of element nodes on both faces
	 * \param field_dim number of field dimensions
	 * \param coords array of nodal coordinates in local ordering.
	 *        The number of nodes in the array is either the number
	 *        of face nodes or the total number of nodes. The class
	 *        handles each case accordingly. */
	SurfaceShapeT(GeometryT::CodeT geometry_code, int num_ip, int num_nodes,
		int num_nodes_per_facet, int field_dim, const LocalArrayT& coords);
		
	/** constructor. 
	 * \param link shared parent domain and "synch-ed" CurrIP and shared
	 *        parent domain.
	 * \param coords array of nodal coordinates in local ordering.
	 *        The number of nodes in the array is either the number
	 *        of face nodes or the total number of nodes. The class
	 *        handles each case accordingly. */	 
	SurfaceShapeT(const SurfaceShapeT& link, const LocalArrayT& coords);

	/** total number of element nodes */
	int TotalNodes(void) const;
	
	/** number of nodes on one face. This is usually just SurfaceShapeT::TotalNodes/2 */
	int NumFacetNodes(void) const;

	/** dimension of the field */
	int FieldDim(void) const;

	/** set all local parameters */
	virtual void Initialize(void);

	/** \name current integration point */
	/*@{*/
	/** interpolate the jump in the field values to the current integration point 
	 * \param nodal array of nodal values: [nen] x [nu]
	 * \return interpolated jump in the nodal values: [nu] */
	const dArrayT& InterpolateJumpU(const LocalArrayT& nodal);

	/** interpolate the jump in the given nodal values to the current integration
	 * point.
	 * \param nodal array of nodal values: [nen] x [nv]
	 * \param jump  interpolated jump in the nodal values: [nv] */
	void InterpolateJump(const LocalArrayT& nodal, dArrayT& jump) const;	

	/** interpolate field values to the current integration point. 
	 * \param nodal array of nodal values. The number of nodal
	 *        values passed in must be either SurfaceShapeT::TotalNodes
	 *        or SurfaceShapeT::NumFacetNodes, that is: [nnd] x [nu].
	 *        If nnd is SurfaceShapeT::NumFacetNodes, then the values
	 *        passed in are assumed to be on the "first" face, i.e.,
	 *        the nodes defined in the first half of the connectivities
	 *        of an element.
	 * \return interpolated jump in the nodal values: [nu] */
	void Interpolate(const LocalArrayT& nodal, dArrayT& u) const;

	/** coordinates of the current integration point */
	const dArrayT& IPCoords(void);

	/** extrapolate integration point values to the nodes.
	 * \param IPvalues values from a single integration point: [numvals]
	 * \param nodalvalues extrapolated values: [fNumNodes] x [numvals] */
	void Extrapolate(const dArrayT& IPvalues, dArray2DT& nodalvalues);

	/** jump gradient table:
	 *
	 *     fgrad_d = d delta_i/d u_j	[i] = FieldDim
	 *                              	[j] = NumNodes*FieldDim */
	const dMatrixT& Grad_d(void) const;

    /** jump gradient table:
	 *
	 *     fgrad_dTgrad_d = d delta_k/d u_i	d delta_k/d u_j
	 *                          	[i],[j] = NumNodes*FieldDim */
	const dMatrixT& Grad_dTGrad_d(void) const;

	/** jacobian of the area transformation using the nodes on the 1st facet */
	double Jacobian(void);

	/** jacobian of the area transformation using the nodes on the 1st facet 
	 * \param Q returns with the transformation from the global coordinates
	 *        to the local coordinates of the element */
	double Jacobian(dMatrixT& Q);

	/** jacobian of the area transformation using the nodes on the 1st facet 
	 * \param Q returns with the transformation from the global coordinates
	 *        to the local coordinates of the element
	 * \param dQ the third rank linearization of the transformation Q with
	 *        respect to coordinates of the surface. Each matrix in the list
	 *        if the linearization of a column vector of Q. */
	double Jacobian(dMatrixT& Q, ArrayT<dMatrixT>& dQ);

	/** nodal shape functions at the current integration point.
	 * \param Na destination for the nodal shape function. Must have length either
	 *        the number of nodes on one face or the number of nodes on both faces. */
	void Shapes(dArrayT& Na) const;

	/** nodal jump shape functions at the current integration point. This is the
	 * this is the number of nodes on both faces */
	void JumpShapes(dArrayT& jump_Na) const { fjumpNa.RowAlias(fCurrIP, jump_Na); };
	/*@}*/

	/** local node numbers on each facet */
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
	int fTotalNodes;   // fNumFacets x fNumFacetNodes
	int fNumFacetNodes;// number nodes on each facet
	int fNumFacets; 
	int fFieldDim;

	/* surface coordinates */
	const LocalArrayT& fCoords;
	LocalArrayT fFacetCoords;

	/** shape functions: [nip] x [n_tot] */
	dArray2DT fNa;	

	/** jump shape functions: [nip] x [n_tot] */
	dArray2DT fjumpNa;	
	
	/* shape function tables */
	ArrayT<dMatrixT> fgrad_d;
	ArrayT<dMatrixT> fgrad_dTgrad_d;
	
	/* shape function derivative tables */
	Array2DT<dMatrixT> fgrad_dd;
	
	/* return value */
	dArrayT fInterp;

	/* return value */
	dArrayT fIPCoord;
	
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

	fDomain->Interpolate(fFacetCoords, fIPCoord, fCurrIP);
	return fIPCoord;
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

inline const dArrayT& SurfaceShapeT::InterpolateJumpU(const LocalArrayT& nodal)
{
	InterpolateJump(nodal, fInterp);
	return fInterp;
}

} // namespace Tahoe 
#endif /* _SURFACE_SHAPE_T_H_ */
