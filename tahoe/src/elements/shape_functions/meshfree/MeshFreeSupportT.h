/* $Id: MeshFreeSupportT.h,v 1.1.1.1 2001-01-29 08:20:33 paklein Exp $ */
/* created: paklein (09/07/1998)                                          */
/* NOTE: Currently, fnNeighborData and feNeighborData don't reflect the   */
/* configuration of the cutting surfaces, so the structure of the         */
/* global equation matrix doesn't need to be reset when the crack         */
/* grows, although new LHS matrices should be computed. In order for      */
/* the global equation matrix to change fnNeighborData and                */
/* would need to be recomputed.                                           */
/* base class for support to MLS shape functions                          */

#ifndef _MF_SUPPORT_T_H_
#define _MF_SUPPORT_T_H_

/* base class */
#include "MeshFreeT.h"

/* direct members */
#include "iArrayT.h"
#include "dArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "RaggedArray2DT.h"
#include "nArrayGroupT.h"
#include "nVariArray2DT.h"

/* forward declarations */
class dArray2DT;
class ParentDomainT;
class OrthoMLSSolverT;
class MLSSolverT;
class iGridManagerT;
class iNodeT;

class MeshFreeSupportT: public MeshFreeT
{
public:

	/* constructor */
	MeshFreeSupportT(const ParentDomainT& domain, const dArray2DT& coords,
		const iArray2DT& connects, const iArrayT& nongridnodes, FormulationT code,
		double dextra, int complete, bool store_shape);

	/* destructor */
	virtual ~MeshFreeSupportT(void);
	
	/* steps to initialization - modifications to the support size must
	 * occur before setting the neighbor data */
	virtual void SetSupportSize(void);
	virtual void SetNeighborData(void);

	/* specify nodes/cells to skip when doing MLS calculations */
	void SetSkipNodes(const iArrayT& skip_nodes);
	const iArrayT& SkipNodes(void) const;
	void SetSkipElements(const iArrayT& skip_elements);
	const iArrayT& SkipElements(void) const;

	/* synchronize Dmax with another set (of active EFG nodes) */
	void SynchronizeDmax(dArrayT& synchDmax);
	void SetDmax(const iArrayT& node, const dArrayT& Dmax);
	void GetDmax(const iArrayT& node, dArrayT& Dmax) const;
	const dArrayT& Dmax(void) const;

	/* cutting facet functions */
	virtual void SetCuttingFacets(const dArray2DT& facet_coords, int num_facet_nodes);
	void ResetFacets(const ArrayT<int>& facets);
	const ArrayT<int>& ResetNodes(void) const; // after passing in new facets, hold
	const ArrayT<int>& ResetCells(void) const; // ids of affected nodes and cells

	/* "load" data for the specified node (global numbering) */
	void LoadNodalData(int node, iArrayT& neighbors, dArrayT& phi,
		dArray2DT& Dphi);

	/* "load" data for the specified element (0...)
	 * for all integration points in the element */
	void LoadElementData(int element, iArrayT& neighbors,
		dArray2DT& phi, ArrayT<dArray2DT>& Dphi);

	/* setting the MLS functions at an arbitrary point - return 1 if successful */
	int SetFieldAt(const dArrayT& x, const dArrayT* shift = NULL);
	int SetFieldUsing(const dArrayT& x, const ArrayT<int>& nodes);
	const ArrayT<int>& NeighborsAt(void) const;
	const dArrayT& FieldAt(void) const;	
	const dArray2DT& DFieldAt(void) const;

	/* collect all nodes covering the point x - returns 0 if not enough found */
	int BuildNeighborhood(const dArrayT& x, AutoArrayT<int>& nodes);

	/* access to neighbors database */
	const iArrayT& ElementNeighborsCounts(void) const;
	const RaggedArray2DT<int>& ElementNeighbors(void) const;
	const RaggedArray2DT<int>& NodeNeighbors(void) const;

	/* list of nodes used in the connectivities */
	const iArrayT& NodesUsed(void) const;

	/* write MLS statistics */
	void WriteStatistics(ostream& out) const;

	/* nodal coordinates */
	const dArray2DT& NodalCoordinates(void) const;

protected:

	/* shape function state */
	enum ShapeState {kNotInit =-1,
	                kNoReform = 0,
	                  kReform = 1};

	/* initialize search grid */
	void SetSearchGrid(void);

	/* generate lists of all nodes that fall within Dmax of the
	 * nodal coords (self included) */
	void SetNodeNeighborData(const dArray2DT& coords);

	/* generate lists of all nodes that fall within Dmax of the
	 * element integration points */
	void SetElementNeighborData(const iArray2DT& connects);

	/* compute all nodal shape functions and derivatives */
	virtual void SetNodalShapeFunctions(void);

	/* compute all integration point shape functions and derivatives */
	virtual void SetElementShapeFunctions(void);

	/* allocate and set pointers for shape function databases */
	virtual void InitNodalShapeData(void);
	virtual void InitElementShapeData(void);

	/* process boundaries - nodes marked as "inactive" at the
	 * current x_node by setting dmax = -1.0 */
	virtual void ProcessBoundaries(const dArray2DT& coords,
		const dArrayT& x_node, dArrayT& dmax) = 0;
	virtual int Visible(const double* x1, const double* x2) = 0;

private:

	/* computing the MLS fits */
	void ComputeElementData(int element, iArrayT& neighbors, dArray2DT& phi,
		ArrayT<dArray2DT>& Dphi);
	void ComputeNodalData(int node, const iArrayT& neighbors, dArrayT& phi,
		dArray2DT& Dphi);

	/* set dmax for each node in the connectivities */
	void SetDmax(void);
	void SetDmaxFromConnects(void); // faster, but not strictly correct
	void SetNodesUsed(void);

	/* swap data */
	void SwapData(const iArrayT& counts, iArray2DT** pfrom, iArray2DT** pto);

protected:

	/* MLS parameters */
	double fDextra;
	int    fComplete;   // order of field completeness
	int    fStoreShape; // 1 => compute and store all shape functions
	                    // 0 => compute shape functions as needed

	/* cutting surface data */
	int fNumFacetNodes;
	const dArray2DT* fCutCoords;

	/* nodal coordinates */
	const dArray2DT& fCoords;

	/* parent integration domain */
	const ParentDomainT& fDomain;

/* MLS solvers */
	OrthoMLSSolverT* fEFG;
	MLSSolverT*      fRKPM;

	/* search grid */
	iGridManagerT* fGrid;

	/* nodal neighbor lists */
	iArrayT fSkipNode;
	iArrayT fnNeighborCount;
	RaggedArray2DT<int> fnNeighborData;
	
	/* element neighbor lists */
	iArrayT fSkipElement;
	iArrayT feNeighborCount;
	RaggedArray2DT<int> feNeighborData;

	/* MLS computation work space */
	AutoArrayT<int>       fneighbors; //used by SetFieldAt
	dArrayT               fvolume;
	dArrayT               fdmax, fdmax_ip;
	nArrayGroupT<double>  fdmax_man;
	dArray2DT             fcoords;
	nVariArray2DT<double> fcoords_man;
	dArray2DT             fx_ip_table;
	dArrayT               felShapespace;
	dArrayT               fndShapespace;

	/* external data */
	const iArray2DT& fConnects; // element connectivities (global numbering)
	const iArrayT&   fNonGridNodes; // EFG nodes not on the integration grid (global numbering)

	/* nodal attributes */
	dArrayT fVolume; // nodal volume (integration weight) -> just 1.0 for now
	dArrayT fDmax;   // nodal support size

	/* nodal neighbor lists */
	iArrayT fNodesUsed; // (ordered) compact list of nodes
	                    // (global numbering) used in the
	                    // EFG domain (connectivities + off-grid nodes)
	
	/* nodal shape function database */
	RaggedArray2DT<double> fnPhiData;
	RaggedArray2DT<double> fnDPhiData;
	
	/* element shape function database */
	RaggedArray2DT<double> fePhiData;
	RaggedArray2DT<double> feDPhiData;
	
	/* selective recompute lists (recomputes all if lists empty) */
	AutoArrayT<int> fResetNodes;
	AutoArrayT<int> fResetElems;

	/* runtime flags */
	ShapeState fReformNode;
	ShapeState fReformElem;

};

/* inlines */
inline const dArrayT& MeshFreeSupportT::Dmax(void) const { return fDmax; }
inline const iArrayT& MeshFreeSupportT::ElementNeighborsCounts(void) const { return feNeighborCount; }
inline const RaggedArray2DT<int>& MeshFreeSupportT::ElementNeighbors(void) const { return feNeighborData; }
inline const RaggedArray2DT<int>& MeshFreeSupportT::NodeNeighbors(void) const { return fnNeighborData; }
inline const iArrayT& MeshFreeSupportT::NodesUsed(void) const { return fNodesUsed; }

inline const iArrayT& MeshFreeSupportT::SkipNodes(void) const { return fSkipNode; }
inline const iArrayT& MeshFreeSupportT::SkipElements(void) const { return fSkipElement; }
inline const dArray2DT& MeshFreeSupportT::NodalCoordinates(void) const { return fCoords; }
inline const ArrayT<int>& MeshFreeSupportT::NeighborsAt(void) const { return fneighbors; }

inline const ArrayT<int>& MeshFreeSupportT::ResetNodes(void) const { return fResetNodes; }
inline const ArrayT<int>& MeshFreeSupportT::ResetCells(void) const { return fResetElems; }

#endif /* _MF_SUPPORT_T_H_ */
