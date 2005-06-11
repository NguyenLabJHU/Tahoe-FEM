/* $Id: CellGeometryT.h,v 1.6 2005-01-28 00:38:02 cjkimme Exp $ */
#ifndef _CELL_GEOMETRY_T_
#define _CELL_GEOMETRY_T_

#include "SCNIMFT.h"
#include "ParameterInterfaceT.h"
#include "ElementSupportT.h"
#include "iArrayT.h"
#include "dArray2DT.h"
#include "RaggedArray2DT.h"
#include "MeshFreeNodalShapeFunctionT.h"
#include "LinkedListT.h"

namespace Tahoe {

class dArrayT;

class CellGeometryT: public ParameterInterfaceT
{
public:

	/** constructor */
	CellGeometryT(const ElementSupportT& support, bool isAxisymmetric);
	CellGeometryT(void);

	/** collect-geometry specific mesh, node, element, or body information */
	virtual void DefineElements(const ArrayT<StringT>& block_ID, const ArrayT<int>& mat_index);
	
	/** set the nodal coordinates and shape functions
	* \param nodes global ID's of the nodes associated with the rows in the coordinate array */
	void SetNodesAndShapes(const iArrayT* nodes, const dArray2DT* nodal_coordinates, 
		MeshFreeNodalShapeFunctionT* nodalShapeFunctions);
	
	/** generate data structures for integration over the body boundary */
	virtual void BoundaryShapeFunctions(RaggedArray2DT<double>& phis, RaggedArray2DT<int>& supports, dArray2DT& normals) = 0;
	
	void SetNodalElements(SCNIMFT* scnimft_group) { fscnimft = scnimft_group; };
	
	/** compute B matrices for strain smoothing/nodal integration */
	virtual void ComputeBMatrices(RaggedArray2DT<int>& nodalCellSupports, RaggedArray2DT<dArrayT>& bVectorArray,
									dArrayT& cellVolumes, dArray2DT& cellCentroids, RaggedArray2DT<double>& circumferential_B) = 0;

	/** accessor to the element support class */
	const ElementSupportT& ElementSupport(void) const;
	
	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;
	
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
	
protected: /* for derived classes only */

	/** Given a list of nodes covering a point, merge that list and its values into
		an accumulated structure that is the union of the new nodes and the existing
		ones in the data structure. This routine is an insert into sorted routine. */	
	void MergeFacetIntegral(int node_num, double weight, dArrayT& facetNormal, const dArrayT& phiValues,
						const iArrayT& neighbors); 
	
	/** Same as merge facet integral, but here only a single value rather than a vector
	    is merged. Also, the insertionQ flag tells whether to actually insert (insertionQ = true)
	    or to overwrite (insertionQ = false). This latter case is relevant to the axisymmetric
	    elements where the shape function values at the node are not computed during loops
	    over facet integration points. */					
	void MergeNodalValues(int node_num, dArrayT& values, const iArrayT& neighbors, 
						ArrayT< LinkedListT<int> >& suppWorkSpace, 
						ArrayT< LinkedListT<double> >& valWorkSpace, bool insertionQ);
	
	/** Move data from linked list workspaces to RaggedArray2DTs. Also finishes computation of circumferential
		components of B-vectors. */
	void ConfigureDataStructures(RaggedArray2DT<int>& cellSupports, RaggedArray2DT<dArrayT>& bVectors,
							RaggedArray2DT<double>& circumferential_B, dArrayT& cellVolumes);
	
	/** number of integration points per facet for cell volume boundary integration */
	int fNumIP; 

	const ElementSupportT* fElementSupport;

	/** global ID of the nodes in CellGeometryT::fNodalCoordinates */
	const iArrayT* fNodes;
	const dArray2DT* fNodalCoordinates;
	
	/** shape functions */
	MeshFreeNodalShapeFunctionT* fNodalShapes;
	
	/** helper monkey class -- has global and local numbering scheme */
	SCNIMFT* fscnimft;
	
	/** axisymmetric? -- false by default */
	bool qIsAxisymmetric;
	
	/** workspaces for computeBMatrices */
	ArrayT< LinkedListT<int> > nodeWorkSpace;
	ArrayT< LinkedListT<dArrayT> > facetWorkSpace;
	ArrayT< LinkedListT<double> > circumferentialWorkSpace;
	
};

/* accessor to the element support class */
inline const ElementSupportT& CellGeometryT::ElementSupport(void) const
{
#if __option(extended_errorcheck)
	if (!fElementSupport) ExceptionT::GeneralFail("CellGeometryT::ElementSupport", "pointer not set");
#endif
	return *fElementSupport;
}

} /* namespace Tahoe */

#endif /* _CELL_GEOMETRY_T */

