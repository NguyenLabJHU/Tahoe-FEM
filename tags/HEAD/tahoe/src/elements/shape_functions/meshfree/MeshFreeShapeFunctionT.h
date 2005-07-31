/* $Id: MeshFreeShapeFunctionT.h,v 1.1.1.1 2001-01-29 08:20:31 paklein Exp $ */
/* created: paklein (09/10/1998)                                          */
/* MLS shape functions for the displacement interpolation - for           */
/* small strain or total Lagrangian finite deformation. DO NOT            */
/* USE for updated Lagrangian finite deformation.                         */

#ifndef _MF_SHAPE_FUNCTION_T_H_
#define _MF_SHAPE_FUNCTION_T_H_

//TEMP
#include <fstream.h>

/* base class */
#include "ShapeFunctionT.h"

/* direct members */
#include "MeshFreeT.h"
#include "iArray2DT.h"
#include "iAutoArrayT.h"
#include "iArrayT.h"

/* forward declarations */
class MeshFreeSupportT;
template <class TYPE> class RaggedArray2DT;

class MeshFreeShapeFunctionT: public ShapeFunctionT
{
public:

/* constructors */
	MeshFreeShapeFunctionT(GeometryT::CodeT geometry_code, int numIP,
		const LocalArrayT& coords, const dArray2DT& all_coords,
		const iArray2DT& connects, const iArrayT& nongridnodes,
		MeshFreeT::FormulationT code, double dextra, int complete, bool store_shape,
		const int& currelement);

	/* destructor */
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

	/* read/write Dmax */
	void SetDmax(const iArrayT& node, const dArrayT& Dmax);
	void GetDmax(const iArrayT& node, dArrayT& Dmax) const;
	const dArrayT& Dmax(void) const;

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

	/* write MLS statistics */
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

	/* MLS database support */
	MeshFreeSupportT* fMFSupport;
	
	/* current element number */
	const int& fCurrElement;
	
	/* ip data loaded from meshfree */
	iArrayT           fNeighbors;
	dArray2DT         fNaU;
	ArrayT<dArray2DT> fDNaU;

	/* interpolant nodes */
	const iArray2DT& fXConnects; // integration grid cell nodes
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
	if (!fMFSupport) throw eGeneralFail;
	return *fMFSupport;
}
inline const iArrayT& MeshFreeShapeFunctionT::Neighbors(void) const { return fNeighbors; }

#endif /* _MF_SHAPE_FUNCTION_T_H_ */
