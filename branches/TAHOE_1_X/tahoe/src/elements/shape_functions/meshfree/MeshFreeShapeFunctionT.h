/* $Id: MeshFreeShapeFunctionT.h,v 1.8 2003-01-27 07:00:30 paklein Exp $ */
/* created: paklein (09/10/1998) */

#ifndef _MF_SHAPE_FUNCTION_T_H_
#define _MF_SHAPE_FUNCTION_T_H_

/* base class */
#include "ShapeFunctionT.h"

/* direct members */
#include "iArray2DT.h"
#include "iAutoArrayT.h"
#include "iArrayT.h"


namespace Tahoe {

/* forward declarations */
class MeshFreeSupportT;
class ifstreamT;
template <class TYPE> class RaggedArray2DT;

/** interface for meshfree shape functions. See ShapeFunctionT
 * for documentation. */
class MeshFreeShapeFunctionT: public ShapeFunctionT
{
public:

	/** constructors
	 * \param geometry_code geometry of the integration cells 
	 * \param numIP number of integration points per cell
	 * \param coords array of cell nodal coordinates in local ordering
	 * \param all_coords reference to the entire coordinate list
	 * \param connectivities of the integration grid and declares on-grid nodes
	 * \param nongridnodes list of nodes not on the grid
	 * \param currelement reference to the current cell of evaluation
	 * \param in input stream */
	MeshFreeShapeFunctionT(GeometryT::CodeT geometry_code, int numIP,
		const LocalArrayT& coords, const dArray2DT& all_coords,
		const iArray2DT& connects, const iArrayT& nongridnodes,
		const int& currelement, ifstreamT& in);

	/** destructor */
	~MeshFreeShapeFunctionT(void);

	/* initialization - modifications to the support size must
	 * occur before setting the neighbor data. Coordinates and
	 * connecitivies must be set */
	void SetSupportSize(void);
	void SetNeighborData(void);
	void SetExactNodes(const iArrayT& exact_nodes);

	/* specify nodes/cells to skip when doing MLS calculations */
	void SetSkipNodes(const iArrayT& skip_nodes);
	void SetSkipElements(const iArrayT& skip_elements);

	/* read/write nodal meshfree parameters */
	void SetNodalParameters(const iArrayT& node, const dArray2DT& nodal_params);
	void GetNodalParameters(const iArrayT& node, dArray2DT& nodal_params) const;
	dArray2DT& NodalParameters(void);

	/* compute global shape derivatives */ 	
	virtual void SetDerivatives(void);
	int SetDerivativesAt(const dArrayT& x); // returns 0 if MLS fails
	void UseDerivatives(const iArrayT& neighbors, const dArray2DT& Dfield); // load external values

	/* cutting facet functions */
	void SetCuttingFacets(const dArray2DT& facet_coords, int num_facet_nodes);
	void ResetFacets(const ArrayT<int>& facets);
	const ArrayT<int>& ResetNodes(void) const;
	const ArrayT<int>& ResetCells(void) const;

	/* returns the number of neighbors for the current element */
	int NumberOfNeighbors(void) const;
	const iArrayT& Neighbors(void) const;

	/* access to MLS field neighbor data */
	const iArrayT& ElementNeighborsCounts(void) const;
	const RaggedArray2DT<int>& ElementNeighbors(void) const;
	const RaggedArray2DT<int>& NodeNeighbors(void) const;

	/* reconstruct displacement field */
	void SelectedNodalField(const dArray2DT& all_DOF, const iArrayT& nodes, dArray2DT& field);
	void NodalField(const dArray2DT& DOF, dArray2DT& field, iArrayT& nodes);
	void NodalField(const dArray2DT& DOF, dArray2DT& field, dArray2DT& Dfield,
		iArrayT& nodes);

	/* print the current ip shape functions to the output stream */
	virtual void Print(ostream& out) const;
	void PrintAt(ostream& out) const;

	/* write MLS information */
	void WriteParameters(ostream& out) const;
	void WriteStatistics(ostream& out) const;

	/* blend FE/MLS shape functions for interpolant nodes */
	void BlendElementData(void);
	void BlendNodalData(int node, const iArrayT& nodes, dArrayT& phi);

	/* reference to the support */
	MeshFreeSupportT& MeshFreeSupport(void) const;

private:

	/* initialize blending database */
	void InitBlend(void);

protected:

	/** MLS database support */
	MeshFreeSupportT* fMFSupport;
	
	/** reference to the current element number */
	const int& fCurrElement;
	
	/* ip data loaded from meshfree */
	iArrayT           fNeighbors;
	dArray2DT         fNaU;
	ArrayT<dArray2DT> fDNaU;

	/** list of integration grid cell connectivities */
	const iArray2DT& fXConnects; 
	iArrayT   fExactNodes; // 1...
	                       // should this be a copy of a reference to
	                       // a dynamically changing list? would have
	                       // to be sure to re-verify new list everytime.
	                       // copy for now.
	iArrayT   fElemHasExactNode;   // flag and map to element data
	iArray2DT fElemFlags;
	
	/* work space for blended shape functions */
	dArrayT   fR;  // blending ramp function
	dArray2DT fDR; // and derivatives
	dArray2DT fNa_tmp;
	ArrayT<dArray2DT> fDNa_tmp;
	dArrayT     felSpace;
	dArrayT     fndSpace;
	iAutoArrayT fNeighExactFlags; // 1 if neighbor is interpolant node
};

/* inlines */
inline MeshFreeSupportT& MeshFreeShapeFunctionT::MeshFreeSupport(void) const
{
	if (!fMFSupport) throw ExceptionT::kGeneralFail;
	return *fMFSupport;
}
inline const iArrayT& MeshFreeShapeFunctionT::Neighbors(void) const { return fNeighbors; }

} // namespace Tahoe 
#endif /* _MF_SHAPE_FUNCTION_T_H_ */
