/* $Id: MeshFreeNodalShapeFunctionT.h,v 1.4 2004-07-29 23:42:12 cjkimme Exp $ */
#ifndef _MF_NODAL_SHAPE_FUNCTION_T_H_
#define _MF_NODAL_SHAPE_FUNCTION_T_H_

/* direct members */
#include "iArray2DT.h"
#include "iAutoArrayT.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "dArray2DT.h"


namespace Tahoe {

/* forward declarations */
class MeshFreeSupportT;
class ParameterListT;
class LocalArrayT;
template <class TYPE> class RaggedArray2DT;

/** interface for meshfree shape functions. See ShapeFunctionT
 * for documentation. */
class MeshFreeNodalShapeFunctionT
{
public:

	/** constructors
	 * \param geometry_code geometry of the integration cells 
	 * \param numIP number of integration points per cell
	 * \param nodes array of cell nodal coordinates in local ordering
	 * \param all_coords reference to the entire coordinate list
	 * \param connectivities of the integration grid and declares on-grid nodes
	 * \param nonNodes evaluate shape functions here, but they're not nodes
	 * \param currelement reference to the current cell of evaluation
	 * \param in input stream */
	MeshFreeNodalShapeFunctionT(int numSD,/* const LocalArrayT& nodes,*/ const dArray2DT& all_coords,
		const iArray2DT& connects, const dArray2DT& nonNodes, const ParameterListT& mf_support_params);

	/** destructor */
	~MeshFreeNodalShapeFunctionT(void);

	/* initialization - modifications to the support size must
	 * occur before setting the neighbor data. Coordinates and
	 * connecitivies must be set */
	void SetSupportSize(void);
	void SetNeighborData(void);

	/* specify nodes/cells to skip when doing MLS calculations */
	void SetSkipNodes(const iArrayT& skip_nodes);

	/* read/write nodal meshfree parameters */
	void SetNodalParameters(const iArrayT& node, const dArray2DT& nodal_params);
	void GetNodalParameters(const iArrayT& node, dArray2DT& nodal_params) const;
	dArray2DT& NodalParameters(void);
	dArrayT& NodalVolumes(void);

	/* compute shape function at arbitrary point */
	virtual int SetFieldAt(const dArrayT& x, const dArrayT* shift); // returns 0 if MLS fails
	const dArrayT& FieldAt(void);

	/* compute global shape derivatives (at arbitrary point)*/ 	
	virtual void SetDerivatives(void);
	int SetDerivativesAt(const dArrayT& x); // returns 0 if MLS fails
	const dArray2DT& DFieldAt(void);
	void UseDerivatives(const iArrayT& neighbors, const dArray2DT& Dfield); // load external values

	/* cutting facet functions */
	void SetCuttingFacets(const dArray2DT& facet_coords, int num_facet_nodes);
	void ResetFacets(const ArrayT<int>& facets);
	const ArrayT<int>& ResetNodes(void) const;
	const ArrayT<int>& ResetCells(void) const;

	/* returns the number of neighbors for the current node */
	int NumberOfNeighbors(void) const;
	const iArrayT& Neighbors(void) const;

	/* access to MLS field neighbor data */
	//const iArrayT& ElementNeighborsCounts(void) const;
	//const RaggedArray2DT<int>& ElementNeighbors(void) const;
	const RaggedArray2DT<int>& NodeNeighbors(void) const;

	/* reconstruct displacement field */
	void SelectedNodalField(const dArray2DT& all_DOF, const iArrayT& nodes, dArray2DT& field);
	void NodalField(const dArray2DT& DOF, dArray2DT& field, iArrayT& nodes);
	void NodalField(const dArray2DT& DOF, dArray2DT& field, dArray2DT& Dfield,
		iArrayT& nodes);

	/* print the current shape functions to the output stream */
	virtual void Print(ostream& out) const;
	void PrintAt(ostream& out) const;

	/* write MLS information */
	void WriteParameters(ostream& out) const;
	void WriteStatistics(ostream& out) const;

	/* reference to the support */
	MeshFreeSupportT& MeshFreeSupport(void) const;

protected:

	/** MLS database support */
	MeshFreeSupportT* fMFSupport;
	
	/** Number of spatial dimensions */
	int fSD;
	
	/** reference to the current node number */
//	const int& fCurrNode;
	
	/* nodal data loaded from meshfree support classes */
	iArrayT           fNeighbors;
	dArray2DT         fNaU;
	ArrayT<dArray2DT> fDNaU;
	
	/** need to hold this to pass a reference to MFSupport classes */
	iArrayT fNonGridNodes;

	/** list of integration grid cell connectivities */
//	const iArray2DT& fXConnects; 
//	iArrayT   fExactNodes; // 1...
	                       // should this be a copy of a reference to
	                       // a dynamically changing list? would have
	                       // to be sure to re-verify new list everytime.
	                       // copy for now.
//	iArrayT   fElemHasExactNode;   // flag and map to element data
//	iArray2DT fElemFlags;
	
};

/* inlines */
inline MeshFreeSupportT& MeshFreeNodalShapeFunctionT::MeshFreeSupport(void) const
{
	if (!fMFSupport) throw ExceptionT::kGeneralFail;
	return *fMFSupport;
}

inline const iArrayT& MeshFreeNodalShapeFunctionT::Neighbors(void) const { return fNeighbors; }

} // namespace Tahoe 
#endif /* _MF_SHAPE_FUNCTION_T_H_ */
