/* $Id: D3MeshFreeSupportT.h,v 1.2 2004-10-30 20:51:19 raregue Exp $ */
/* created: paklein (10/23/1999) */
#ifndef _D3_MF_SUPPORT_T_H_
#define _D3_MF_SUPPORT_T_H_

/* base class */
#include "D2MeshFreeSupportT.h"

namespace Tahoe {

/** class for support of meshfree field calculations up to
 * third gradients. See documentation for D2MeshFreeSupportT 
 * for information about class initialization. */
class D3MeshFreeSupportT: public D2MeshFreeSupportT
{
public:

	/** constructor.
	 * \param domain used to determine the location of integration points
	 * \param coords array of all particle coordinates 
	 * \param connects integration cell connectivities 
	 * \param nongridnodes index of paricles not included in the connectivities
	 * \param in input stream for class and window function parameters */
	D3MeshFreeSupportT(const ParentDomainT* domain, const dArray2DT& coords,
		const iArray2DT& connects, const iArrayT& nongridnodes);

	/** determine nodal support parameters based window function parameters */
	virtual void InitNeighborData(void);

	/** fetch data for the specified node. triggers recalculation of any
	 * nodal shape function that have been reset.
	 * \param node particle data to fetch
	 * \param neighbors returns with neighbors of node: [nnd]
	 * \param phi returns with the shape function values of neighbors at node: [nnd]
	 * \param Dphi returns with neighbors shape function derivatives at node: [nsd] x [nnd]
	 * \param DDphi returns with neighbors shape function second derivatives at node: [nstr] x [nnd]
	 * \param DDDphi returns with neighbors shape function second derivatives at node: [nsd*nsd] x [nnd] */	
	void LoadNodalData(int node, iArrayT& neighbors, dArrayT& phi,
		dArray2DT& Dphi, dArray2DT& DDphi, dArray2DT& DDDphi);

	/** fetch data for the specified integration cell. triggers recalculation of any
	 * integration cell shape functions that have been reset.
	 * \param element cell data to fetch
	 * \param neighbors returns with neighbors of all integration points in the cell: [nnd]
	 * \param phi returns with the shape function values of neighbors at the integration points: [nip] x [nnd]
	 * \param Dphi returns with neighbor shape function derivatives: [nip] x [nsd] x [nnd]
	 * \param DDphi returns with neighbor shape function second derivatives: [nip] x [nstr] x [nnd]
	 * \param DDDphi returns with neighbor shape function third derivatives: [nip] x [nsd*nsd] x [nnd] */
	void LoadElementData(int element, iArrayT& neighbors,
		dArray2DT& phi, ArrayT<dArray2DT>& Dphi, ArrayT<dArray2DT>& DDphi, ArrayT<dArray2DT>& DDDphi);

	/** set field at x.
	 * \param x arbitrary field point
	 * \param shift pointer to shift of x to use when determining the neighbor particles. NULL
	          for no shift.
     * \return 1 if successful, 0 otherwise */
	//inherited
	//int SetFieldAt(const dArrayT& x, const dArrayT* shift = NULL);

	/** shape function third derivatives for NeighborsAt the last call to 
	 * SetFieldAt or SetFieldUsing
	 * \return 2D array dimension: [nsd*nsd] x [nnd] */
	const dArray2DT& DDDFieldAt(void) const;

protected:

	/* compute all nodal shape functions and derivatives */
	virtual void SetNodalShapeFunctions(void);

	/* compute all integration point shape functions and derivatives */
	virtual void SetElementShapeFunctions(void);

	/* allocate and set pointers for shape function databases */
	void InitNodalShapeData(void);
	void InitElementShapeData(void);

private:

	/* computing the MLS fits */
	void ComputeNodalData(int node, const iArrayT& neighbors, dArrayT& phi,
		dArray2DT& Dphi, dArray2DT& DDphi, dArray2DT& DDDphi);

	void ComputeElementData(int element, iArrayT& neighbors, dArray2DT& phi,
		ArrayT<dArray2DT>& Dphi, ArrayT<dArray2DT>& DDphi, ArrayT<dArray2DT>& DDDphi);

private:
	
	/* nodal shape function database */
	RaggedArray2DT<double> fnDDDPhiData;
	
	/* element shape function database */
	RaggedArray2DT<double> feDDDPhiData;	
};

} // namespace Tahoe 
#endif /* _D3_MF_SUPPORT_T_H_ */
