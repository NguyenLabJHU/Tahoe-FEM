/* $Id: D2MeshFreeSupportT.h,v 1.10 2004-12-24 20:32:45 kyonten Exp $ */
/* created: paklein (10/23/1999) */
#ifndef _D2_MF_SUPPORT_T_H_
#define _D2_MF_SUPPORT_T_H_

/* base class */
#include "MeshFreeSupportT.h"

namespace Tahoe {

/* forward declarations */
class D2OrthoMLSSolverT;

/** class for support of meshfree field calculations up to
 * second gradients. See documentation for MeshFreeSupportT 
 * for information about class initialization. */
class D2MeshFreeSupportT: public MeshFreeSupportT
{
public:

	/** constructor.
	 * \param domain used to determine the location of integration points
	 * \param coords array of all particle coordinates 
	 * \param connects integration cell connectivities 
	 * \param nongridnodes index of paricles not included in the connectivities
	 * \param in input stream for class and window function parameters */
	D2MeshFreeSupportT(const ParentDomainT* domain, const dArray2DT& coords,
		const iArray2DT& connects, const iArrayT& nongridnodes);

	//************************************************************************************//
	// kyonten
	/** construct object sufficient for calling methods inherited from ParameterInterfaceT
	 * to collect the class parameters, but not for doing any meshfree calculations */
	D2MeshFreeSupportT(void);
	//************************************************************************************//
	
	/** determine nodal support parameters based window function parameters */
	virtual void InitNeighborData(void);

	/** fetch data for the specified node. triggers recalculation of any
	 * nodal shape function that have been reset.
	 * \param node particle data to fetch
	 * \param neighbors returns with neighbors of node: [nnd]
	 * \param phi returns with the shape function values of neighbors at node: [nnd]
	 * \param Dphi returns with neighbors shape function derivatives at node: [nsd] x [nnd]
	 * \param DDphi returns with neighbors shape function second derivatives at node: [nstr] x [nnd] */
	void LoadNodalData(int node, iArrayT& neighbors, dArrayT& phi,
		dArray2DT& Dphi, dArray2DT& DDphi);

	/** fetch data for the specified integration cell. triggers recalculation of any
	 * integration cell shape functions that have been reset.
	 * \param element cell data to fetch
	 * \param neighbors returns with neighbors of all integration points in the cell: [nnd]
	 * \param phi returns with the shape function values of neighbors at the integration points: [nip] x [nnd]
	 * \param Dphi returns with neighbor shape function derivatives: [nip] x [nsd] x [nnd]
	 * \param DDphi returns with neighbor shape function second derivatives: [nip] x [nstr] x [nnd] */
	void LoadElementData(int element, iArrayT& neighbors,
		dArray2DT& phi, ArrayT<dArray2DT>& Dphi, ArrayT<dArray2DT>& DDphi);

	/** set field at x.
	 * \param x arbitrary field point
	 * \param shift pointer to shift of x to use when determining the neighbor particles. NULL
	          for no shift.
     * \return 1 if successful, 0 otherwise */
	int SetFieldAt(const dArrayT& x, const dArrayT* shift = NULL);

	/** shape function second derivatives for NeighborsAt the last call to 
	 * SetFieldAt or SetFieldUsing
	 * \return 2D array dimension: [nstr] x [nnd] */
	const dArray2DT& DDFieldAt(void) const;
	
	//*****************************************************************//
	// kyonten: to inherit the ParameterInterfaceT from MeshFreeSupportT
	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
	//*****************************************************************//

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

protected:

	/* meshfree type */
	//FormulationT fMeshfreeType; //kyonten
	
	/* higher order MLS solver */
	D2OrthoMLSSolverT* fD2EFG;
	
	/* nodal shape function database */
	RaggedArray2DT<double> fnDDPhiData;
	
	/* element shape function database */
	RaggedArray2DT<double> feDDPhiData;	
};

} // namespace Tahoe 
#endif /* _D2_MF_SUPPORT_T_H_ */
