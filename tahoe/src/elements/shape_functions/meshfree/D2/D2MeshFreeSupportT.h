/* $Id: D2MeshFreeSupportT.h,v 1.1.1.1 2001-01-29 08:20:33 paklein Exp $ */
/* created: paklein (10/23/1999)                                          */

#ifndef _D2_MF_SUPPORT_T_H_
#define _D2_MF_SUPPORT_T_H_

/* base class */
#include "MeshFreeSupportT.h"

/* forward declarations */
class D2OrthoMLSSolverT;

class D2MeshFreeSupportT: public MeshFreeSupportT
{
public:

	/* constructor */
	D2MeshFreeSupportT(const ParentDomainT& domain, const dArray2DT& coords,
		const iArray2DT& connects, const iArrayT& nongridnodes,  FormulationT code,
		double dextra, int complete, bool store_shape);

	/* steps to initialization - modifications to the support size must
	 * occur before setting the neighbor data */
	virtual void SetNeighborData(void);

	/* "load" data for the specified node (global numbering) */
	void LoadNodalData(int node, iArrayT& neighbors, dArrayT& phi,
		dArray2DT& Dphi, dArray2DT& DDphi);

	/* "load" data for the specified element (0...)
	 * for all integration points in the element */
	void LoadElementData(int element, iArrayT& neighbors,
		dArray2DT& phi, ArrayT<dArray2DT>& Dphi, ArrayT<dArray2DT>& DDphi);

	/* setting the MLS functions at an arbitrary point x */
	int SetFieldAt(const dArrayT& x, AutoArrayT<int>& nodes); // returns 0 if MLS fit fails
		//TEMP - needs to be virtual???

	const dArray2DT& DDFieldAt(void) const;

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
		dArray2DT& Dphi, dArray2DT& DDphi);

	void ComputeElementData(int element, iArrayT& neighbors, dArray2DT& phi,
		ArrayT<dArray2DT>& Dphi, ArrayT<dArray2DT>& DDphi);

private:

/* higher order MLS solver */
	D2OrthoMLSSolverT* fD2EFG;
	
	/* nodal shape function database */
	RaggedArray2DT<double> fnDDPhiData;
	
	/* element shape function database */
	RaggedArray2DT<double> feDDPhiData;	
};

#endif /* _D2_MF_SUPPORT_T_H_ */
