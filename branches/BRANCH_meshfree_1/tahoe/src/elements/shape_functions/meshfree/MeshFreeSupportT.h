/* $Id: MeshFreeSupportT.h,v 1.1.1.1.4.2 2001-06-19 08:58:44 paklein Exp $ */
/* created: paklein (09/07/1998)                                          */

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
#include "nArray2DGroupT.h"
#include "VariArrayT.h"
#include "nVariArray2DT.h"
#include "WindowT.h"

/* forward declarations */
class ifstreamT;
class ofstreamT;
class dArray2DT;
class ParentDomainT;
class OrthoMLSSolverT;
class MLSSolverT;
class iGridManagerT;
class iNodeT;

/** Base class for support of meshfree methods.
 * This class sits between the shape function and the MLS solver. This
 * class feeds the MLS solver coordinates and meshfree parameters for
 * a given field point and can store the resulting shape functions and
 * derivative. Shape functions and there derivatives can subsequently
 * be retrieved from storage or calculated as needed.
 *
 * \note Currently, fnNeighborData and feNeighborData don't reflect the
 * configuration of the cutting surfaces, so the structure of the
 * global equation matrix doesn't need to be reset when the crack
 * grows, although new LHS matrices should be computed. In order for
 * the global equation matrix to change fnNeighborData and would need to be 
 * recomputed. */
class MeshFreeSupportT: public MeshFreeT
{
public:

	/* constructor */
	MeshFreeSupportT(const ParentDomainT& domain, const dArray2DT& coords,
		const iArray2DT& connects, const iArrayT& nongridnodes, ifstreamT& in);

	/* destructor */
	virtual ~MeshFreeSupportT(void);
	
	/* write parameters */
	virtual void WriteParameters(ostream& out) const;
	
	/* steps to initialization - modifications to the support size must
	 * occur before setting the neighbor data */
	virtual void SetSupportSize(void);
	virtual void SetNeighborData(void);

	/* specify nodes/cells to skip when doing MLS calculations */
	void SetSkipNodes(const iArrayT& skip_nodes);
	const iArrayT& SkipNodes(void) const;
	void SetSkipElements(const iArrayT& skip_elements);
	const iArrayT& SkipElements(void) const;

	/* read/write nodal meshfree parameters */
	void SynchronizeNodalParameters(dArray2DT& nodal_params);
	void SetNodalParameters(const iArrayT& node, const dArray2DT& nodal_params);
	void GetNodalParameters(const iArrayT& node, dArray2DT& nodal_params) const;
	const dArray2DT& NodalParameters(void) const;

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
	 * current x_node by setting nodal_params = -1.0 */
	virtual void ProcessBoundaries(const dArray2DT& coords,
		const dArrayT& x_node, dArray2DT& nodal_params) = 0;
	virtual int Visible(const double* x1, const double* x2) = 0;

private:

	/* computing the MLS fits */
	void ComputeElementData(int element, iArrayT& neighbors, dArray2DT& phi,
		ArrayT<dArray2DT>& Dphi);
	void ComputeNodalData(int node, const iArrayT& neighbors, dArrayT& phi,
		dArray2DT& Dphi);

	/* compute nodal support parameters */
	void SetSupport_Spherical_Search(void);
	void SetSupport_Spherical_Connectivities(void); // faster, but not strictly correct
	void SetSupport_Cartesian_Connectivities(void); // faster, but not strictly correct
	void SetNodesUsed(void);

	/* swap data */
	void SwapData(const iArrayT& counts, iArray2DT** pfrom, iArray2DT** pto);

protected:

	/* common meshfree parameters */
	double       fDextra; //used by EFG only
	bool         fStoreShape;
	FormulationT fMeshfreeType;

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
	VariArrayT<double>    fvolume_man;
	dArray2DT             fnodal_param, fnodal_param_ip;	
	nArray2DGroupT<double> fnodal_param_man;
	dArray2DT             fcoords;
	nVariArray2DT<double> fcoords_man;
	dArray2DT             fx_ip_table;
	dArrayT               felShapespace;
	dArrayT               fndShapespace;

	/* external data */
	const iArray2DT& fConnects; // element connectivities (global numbering)
	const iArrayT&   fNonGridNodes; // EFG nodes not on the integration grid (global numbering)

	/* nodal attributes */
	dArrayT fVolume;            // nodal volume (integration weight) -> just 1.0 for now
	dArray2DT fNodalParameters; // nodal meshfree parameters, i.e., the support size

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
inline const dArray2DT& MeshFreeSupportT::NodalParameters(void) const { return fNodalParameters; }
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
