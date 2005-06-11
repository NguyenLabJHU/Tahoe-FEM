/* $Id: MeshFreeSupportT.h,v 1.14 2004-07-22 08:24:11 paklein Exp $ */
/* created: paklein (09/07/1998) */
#ifndef _MF_SUPPORT_T_H_
#define _MF_SUPPORT_T_H_

/* base classes */
#include "MeshFreeT.h"
#include "ParameterInterfaceT.h"

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

namespace Tahoe {

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
 * derivative. Shape functions and their derivatives can subsequently
 * be retrieved from storage or calculated as needed.
 *
 * Initialization of the class involves 2 steps:
 * (I) setting the support parameters for every node:
 *     (1) InitSupportParameters - sets support based on window function requirements
 *     (2) SetSupportParameters - overrides any previously defined values
 *     (3) SynchronizeSupportParameters - take "best" support parameters
 * (II) setting nodal neighbor data - this is only needed if shape function data
 *      is going to be stored.
 *
 * \note Currently, fnNeighborData and feNeighborData don't reflect the
 * configuration of the cutting surfaces, so the structure of the
 * global equation matrix doesn't need to be reset when the crack
 * grows, although new LHS matrices should be computed. In order for
 * the global equation matrix to change fnNeighborData and would need to be 
 * recomputed. */
class MeshFreeSupportT: public MeshFreeT, public ParameterInterfaceT
{
public:

	/** constructor.
	 * \param domain used to determine the location of integration points
	 * \param coords array of all particle coordinates 
	 * \param connects integration cell connectivities 
	 * \param nongridnodes index of paricles not included in the connectivities
	 * \param in input stream for class and window function parameters */
	MeshFreeSupportT(const ParentDomainT* domain, const dArray2DT& coords,
		const iArray2DT& connects, const iArrayT& nongridnodes);

	/** construct object sufficient for calling methods inherited from ParameterInterfaceT
	 * to collect the class parameters, but not for doing any meshfree calculations */
	MeshFreeSupportT(void);

	/** destructor */
	virtual ~MeshFreeSupportT(void);
	
	/** write parameters to out */
	virtual void WriteParameters(ostream& out) const;
	
	/** determine nodal support parameters based window function parameters */
	virtual void InitSupportParameters(void);

	/** set neighbor lists for all field points (particles and integration points) */
	virtual void InitNeighborData(void);

	/** set nodes at which the field should not be calculated */	
	void SetSkipNodes(const iArrayT& skip_nodes); 
	const iArrayT& SkipNodes(void) const; /**< list of nodes to skip */ 

	/** cells whose integration points should not be treated as field points */
	void SetSkipElements(const iArrayT& skip_elements);
	const iArrayT& SkipElements(void) const; /**< list of cells to skip */ 

	/** set support parameters using external data. modify both the internal
	 * support parameters and nodal_params so that both hold the "best" of each. 
	 * \param nodal_params external list of support parameters: [nnd] x [nparam] */
	void SynchronizeSupportParameters(dArray2DT& nodal_params);
	
	/** overwrite nodal support parameters.
	 * \param node list of particles to overwrite : [n] 
	 * \param nodal_params external support parameters: [n] x [nparam] */
	void SetSupportParameters(const iArrayT& node, const dArray2DT& nodal_params);

	/** collent nodal support parameters.
	 * \param node list of particles to collect : [n] 
	 * \param nodal_params support parameters: [n] x [nparam] */
	void GetSupportParameters(const iArrayT& node, dArray2DT& nodal_params) const;
	
	/** read/write access to the nodal parameters */
	dArray2DT& NodalParameters(void);
	
	/** read/write access to nodal integration weights */
	dArrayT& NodalVolumes(void);

	/** set field cutting facets. 
	 * \param facet_coords list of coordinate for each facet: [nfacets] x [num_facet_nodes*nsd] 
	 * \param num_facet_nodes number of nodes defining each facet */
	virtual void SetCuttingFacets(const dArray2DT& facet_coords, int num_facet_nodes);

	/** local stored shape function recalculation.
	 * \param facets facets for which the stored field values should be recalculated */
	void ResetFacets(const ArrayT<int>& facets);
	const ArrayT<int>& ResetNodes(void) const; /**< list of nodes affected by the last call to ResetFacets */
	const ArrayT<int>& ResetCells(void) const; /**< list of cells affected by the last call to ResetFacets */

	/** \name retrieving stored shape function values */
	/*@{*/
	/** fetch data for the specified node. triggers recalculation of any
	 * nodal shape function that have been reset.
	 * \param node particle data to fetch
	 * \param neighbors returns with neighbors of node: [nnd]
	 * \param phi returns with the shape function values of neighbors at node: [nnd]
	 * \param Dphi returns with neighbors shape function derivatives at node: [nsd] x [nnd] */
	void LoadNodalData(int node, iArrayT& neighbors, dArrayT& phi, dArray2DT& Dphi);

	/** fetch data for the specified integration cell. triggers recalculation of any
	 * integration cell shape functions that have been reset.
	 * \param element cell data to fetch
	 * \param neighbors returns with neighbors of all integration points in the cell: [nnd]
	 * \param phi returns with the shape function values of neighbors at the integration points: [nip] x [nnd]
	 * \param Dphi returns with neighbor shape function derivatives: [nip] x [nsd] x [nnd] */
	void LoadElementData(int element, iArrayT& neighbors, dArray2DT& phi, ArrayT<dArray2DT>& Dphi);
	/*@}*/

	/** \name setting the MLS functions at an arbitrary point */
	/*@{*/
	/** set field at x.
	 * \param x arbitrary field point
	 * \param shift pointer to shift of x to use when determining the neighbor particles. NULL
	          for no shift.
     * \return 1 if successful, 0 otherwise */
	int SetFieldAt(const dArrayT& x, const dArrayT* shift = NULL);

	/** set field at x using the specified neighboring particles.
	 * \param x arbitrary field point
	 * \param nodes list of neighboring particles.
     * \return 1 if successful, 0 otherwise */
	int SetFieldUsing(const dArrayT& x, const ArrayT<int>& nodes);

	/** neighbors from the last call to SetFieldAt or SetFieldUsing */
	const ArrayT<int>& NeighborsAt(void) const;

	/** shape function values for NeighborsAt the last call to SetFieldAt or SetFieldUsing
	 * \return array length: [nnd] */
	const dArrayT& FieldAt(void) const;	

	/** shape function derivatives for NeighborsAt the last call to SetFieldAt or SetFieldUsing
	 * \return 2D array dimension: [nsd] x [nnd] */
	const dArray2DT& DFieldAt(void) const;
	/*@}*/

	/** collect all nodes covering an arbitrary point 
	 * \param x arbitrary field point
	 * \param nodes returns with the list of covering nodes
	 * \returns 1 if successful, 0 otherwise */
	int BuildNeighborhood(const dArrayT& x, AutoArrayT<int>& nodes);

	/* access to neighbors database */
	const iArrayT& ElementNeighborsCounts(void) const; /**< list of total neighbors per cell */
	const RaggedArray2DT<int>& ElementNeighbors(void) const; /**< cell neighborhoods */
	const RaggedArray2DT<int>& NodeNeighbors(void) const;    /**< nodal neighborhoods */

	/** list of nodes used in the connectivities */
	const iArrayT& NodesUsed(void) const;

	/** write MLS statistics */
	void WriteStatistics(ostream& out) const;

	/** nodal coordinates */
	const dArray2DT& NodalCoordinates(void) const;

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

	/** construct a new MLSSolverT with the given parameter list. The
	 * MLSSolverT is constructed but MLSSolverT::Initialize still needs to
	 * be called. Host is responsible for freeing the object that is returned. */
	static MLSSolverT* New_MLSSolverT(int nsd, int completeness, const ParameterListT& params);

protected:

	/** state of shape function database */
	enum ShapeState {kNotInit =-1, /**< not initialized */
	                kNoReform = 0, /**< no recalculation needed */
	                  kReform = 1  /**< recalculation needed */ };

	/** initialize search grid */
	void SetSearchGrid(void);

	/* generate lists of all nodes that fall within Dmax of the
	 * nodal coords (self included) */
	void SetNodeNeighborData(const dArray2DT& coords);
	void SetNodeNeighborData_2(const dArray2DT& coords);

	/* generate lists of all nodes that fall within Dmax of the
	 * element integration points */
	void SetElementNeighborData(const iArray2DT& connects);
	void SetElementNeighborData_2(const iArray2DT& connects);

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

	/* test coverage - temporary function until OrthoMLSSolverT class
	 * is brought up to date */
	bool Covers(const dArrayT& field_x, const dArrayT& node_x, int node) const;

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
	const dArray2DT* fCoords;

	/* parent integration domain and its data */
	const ParentDomainT* fDomain;
	int fSD;
	int fIP;

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

	/* pointers to external data */
	const iArray2DT* fConnects; // element connectivities (global numbering)
	const iArrayT*   fNonGridNodes; // EFG nodes not on the integration grid (global numbering)

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
inline dArray2DT& MeshFreeSupportT::NodalParameters(void) { return fNodalParameters; }
inline dArrayT& MeshFreeSupportT::NodalVolumes(void) { return fVolume; }
inline const iArrayT& MeshFreeSupportT::ElementNeighborsCounts(void) const { return feNeighborCount; }
inline const RaggedArray2DT<int>& MeshFreeSupportT::ElementNeighbors(void) const { return feNeighborData; }
inline const RaggedArray2DT<int>& MeshFreeSupportT::NodeNeighbors(void) const { return fnNeighborData; }
inline const iArrayT& MeshFreeSupportT::NodesUsed(void) const { return fNodesUsed; }

inline const iArrayT& MeshFreeSupportT::SkipNodes(void) const { return fSkipNode; }
inline const iArrayT& MeshFreeSupportT::SkipElements(void) const { return fSkipElement; }
inline const dArray2DT& MeshFreeSupportT::NodalCoordinates(void) const { 
	if (!fCoords) ExceptionT::GeneralFail("MeshFreeSupportT::NodalCoordinates", "coordinates not set");
	return *fCoords; 
}
inline const ArrayT<int>& MeshFreeSupportT::NeighborsAt(void) const { return fneighbors; }

inline const ArrayT<int>& MeshFreeSupportT::ResetNodes(void) const { return fResetNodes; }
inline const ArrayT<int>& MeshFreeSupportT::ResetCells(void) const { return fResetElems; }

} // namespace Tahoe 
#endif /* _MF_SUPPORT_T_H_ */
