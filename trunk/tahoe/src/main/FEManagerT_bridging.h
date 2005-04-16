/* $Id: FEManagerT_bridging.h,v 1.25 2005-04-16 02:03:26 paklein Exp $ */
#ifndef _FE_MANAGER_BRIDGING_H_
#define _FE_MANAGER_BRIDGING_H_

/* element configuration header */
#include "ElementsConfig.h"
#ifdef BRIDGING_ELEMENT

/* base class */
#include "FEManagerT.h"

/* direct members */
#include "PointInCellDataT.h"
#include "nMatrixT.h"
#include "KBC_CardT.h"
#include "dArrayT.h"
#include "ElementCardT.h"

namespace Tahoe {

/* forward declarations */
class ParticlePairT;
class BridgingScaleT;
class KBC_ControllerT;
class dSPMatrixT;
class EAMFCC3D;
class EAMT;
class LAdMatrixT; //TEMP

/** extension of FEManagerT for bridging scale calculations */
class FEManagerT_bridging: public FEManagerT
{
public:

	/** constructor */
	FEManagerT_bridging(const StringT& input, ofstreamT& output, CommunicatorT& comm,
		const ArrayT<StringT>& argv, TaskT task);

	/** destructor */
	~FEManagerT_bridging(void);

	/** \name solution update */
	/*@{*/
	/** compute RHS-side, residual force vector and assemble to solver. Aside from
	 * calling FEManagerT::FormRHS, this call also assembles any contributions set
	 * with FEManagerT_bridging::SetExternalForce.
	 * \param group equation group to solve */
	virtual void FormRHS(int group) const;

	/** send update of the solution to the NodeManagerT */
	virtual void Update(int group, const dArrayT& update);

	/** reset the cumulative update vector */
	void ResetCumulativeUpdate(int group);

	/** return the cumulative update. The total update since the last call to
	 * FEManagerT_bridging::ResetCumulativeUpdate */
	const dArrayT& CumulativeUpdate(int group) const { return fCumulativeUpdate[group]; };

	/** set pointer to an external force vector or pass NULL to clear. The array
	 * the length of the number of unknowns for the given group. */
	void SetExternalForce(int group, const dArrayT& external_force);
	
	//DEF added this because I couldn't get the null pointer passed in right - clears the coupling forces
	void ClearExternalForce(int group);
	
	//DEF added to return the pointers to the external force array
	const dArrayT* GetExternalForce(int group) const { return fExternalForce[group]; };
	const dArray2DT* GetExternalForce2D(int group) const { return fExternalForce2D[group]; };

	/** set pointer to an external force vector for the given field */
	void SetExternalForce(const StringT& field, const dArray2DT& external_force, const iArrayT& activefenodes);
	/*@}*/

	/** \name ghost nodes 
	 * The ghost node database must be initialized by calling
	 * FEManagerT_bridging::InitGhostNodes before accessing the lists.*/
	/*@{*/
	/** initialize the ghost node information. Ghost nodes are nodes with externally
	 * prescribed motion. The motion of the ghost nodes must be defined by calling
	 * FEManagerT_bridging::SetGhostNodeKBC.
	 * \param include_image_nodes flag to indicate whether image nodes should be
	 *        included in the list of non-ghost nodes */
	void InitGhostNodes(const StringT& field, const ArrayT<StringT>& ghost_id_list, bool include_image_nodes);

	/** prescribe the motion of ghost nodes. Generate KBC cards to control the
	 * ghost node motion. Assumes all components of the ghost node motion are
	 * prescribed, and that all are prescribed with the same KBC_CardT::CodeT. 
	 * Ghost node information must be initialized by calling 
	 * FEManagerT_bridging::InitGhostNodes first. */
	void SetGhostNodeKBC(KBC_CardT::CodeT code, const dArray2DT& values);

	/** return list of ghost nodes */
	const iArrayT& GhostNodes(void) const { return fGhostNodes; };

	/** return list of ghost nodes */
	const iArrayT& NonGhostNodes(void) const { return fNonGhostNodes; };

	/** compute the ghost-nonghost part of the stiffness matrix */
	void Form_G_NG_Stiffness(const StringT& field, int element_group, dSPMatrixT& K_G_NG);
	/*@}*/

	/** write field values for the given nodes */
	void SetFieldValues(const StringT& field, const iArrayT& nodes, int order,
		const dArray2DT& values);

	/** \name interpolation and projection operators */
	/*@{*/
	/** return the "lumped" (scalar) mass associated with the given nodes */
	void LumpedMass(const iArrayT& nodes, dArrayT& mass) const;
	
	/** initialize interpolation data. Initialize data structures needed to interpolate
	 * field values to the given list of points. Requires that this FEManagerT has
	 * a BridgingScaleT in its element list. */
	void InitInterpolation(const StringT& field, const iArrayT& nodes,
		const dArray2DT& coordinates);

	/** field interpolations. Interpolate the field to the nodes initialized
	 * with the latest call to FEManagerT_bridging::InitInterpolation. */
	void InterpolateField(const StringT& field, int order, dArray2DT& nodal_values);

	/** return the interpolation matrix associated with the active degrees
	 * of freedom */
	void InterpolationMatrix(const StringT& field, dSPMatrixT& G_Interpolation) const;

	/** compute global interpolation matrix for all nodes whose support intersects the MD 
	 *  region, i.e. N_{I}(X_{\alpha}) */
	void Ntf(dSPMatrixT& ntf, const iArrayT& atoms, iArrayT& activefenodes) const;

	/** compute the product with transpose of the interpolation matrix 
	 * \param N data for interpolating from rows in NTf to rows in f
	 * \param f multi-vector in columns
	 * \param f_rows rows if f to include in the product
	 * \param NTf returns with matrix-vector product added to the vectors in columns. Global ids
	 *        in N map directly to rows in NTf.
	 */
	void MultNTf(const PointInCellDataT& N, const dArray2DT& f, const iArrayT& f_rows, dArray2DT& NTf) const;

	/** compute the product with transpose of the interpolation matrix 
	 * \param N data for interpolating from rows in NTf to rows in f
	 * \param f multi-vector in columns
	 * \param f_rows rows if f to include in the product
	 * \param NTf returns with matrix-vector product added to the vectors in columns. Global ids
	 *        in N map directly to rows in NTf.
	 */
	void MultNTf(const InterpolationDataT& N, const dArray2DT& f, const iArrayT& f_rows, dArray2DT& NTf) const;

	/** initialize projection data. Initialize data structures needed to project
	 * field values to the given list of points. Requires that this FEManagerT has
	 * a BridgingScaleT in its element list. */
	void InitProjection(const StringT& field, CommManagerT& comm, const iArrayT& nodes,
		NodeManagerT& node_manager, bool make_inactive);

	/** indicate whether image nodes should be included in the projection */
	virtual bool ProjectImagePoints(void) const;

	/** projects non-field point values onto the mesh. 
	 */
	void Project(const dArray2DT& fine_values, dArray2DT& nodal_values);

	/** project the point values onto the mesh. Project to the nodes using
	 * projection initialized with the latest call to FEManagerT_bridging::InitProjection. */
	void ProjectField(const StringT& field, const NodeManagerT& node_manager, int order);

	/** compute the coarse scale projection at the source points. Project the solution to the source
	 * points initialized with the latest call to FEManagerT_bridging::InitProjection. In other words,
	 * filter out the fine scale part of the solution. */
	void CoarseField(const StringT& field, const NodeManagerT& node_manager, int order, dArray2DT& coarse);

	/** calculate the fine scale part of MD solution as well as total displacement u.  Does not
	  * write into the displacement field */
	void BridgingFields(const StringT& field, NodeManagerT& atom_node_manager,
		NodeManagerT& fem_node_manager, dArray2DT& totalu);
	
	/** calculate the initial FEM displacement via projection of initial MD displacement.  Differs 
	  * from BridgingFields in that projected FE nodal values written into displacement field */
	void InitialProject(const StringT& field, NodeManagerT& atom_node_manager, dArray2DT& projectedu,
		int order);

	/** transpose follower cell data */
	void TransposeFollowerCellData(InterpolationDataT& transpose);

	/** interpolation data */
	const PointInCellDataT& InterpolationData(void) { return fFollowerCellData;	};
	
	/** projection data */
	const PointInCellDataT& ProjectionData(void) { return fDrivenCellData; };
	
	/** list of nodes with projected solutions */
	const iArrayT& ProjectedNodes(void) const { return fProjectedNodes; };
	/*@}*/
	
	/** \name bond density corrections */
	/*@{*/
	/** compute internal correction for the overlap region accounting for ghost node bonds
	 * extending into the coarse scale region */
	void CorrectOverlap_1(const RaggedArray2DT<int>& neighbors, const dArray2DT& coords, double smoothing, double k2);

	/** solve bond densities one at a time enforcing constraints on bond densities using augmented
	 * Lagrange multipliers */
	void CorrectOverlap_2(const RaggedArray2DT<int>& neighbors, const dArray2DT& coords, 
		const StringT& overlap_file, double smoothing, double k2, double k_r, double bound_0, int nip);

	/** solve bond densities one at a time enforcing constraints on bond densities using a
	 * penalty method */
	void CorrectOverlap_22(const RaggedArray2DT<int>& neighbors, const dArray2DT& coords,
		const StringT& overlap_file, double smoothing, double k2, double bound_tol, double stiffness_jump, int nip);

	/** solve bond densities one at a time activating unknowns at integration points only if
	 * the associated domain contains the terminus of a ghost node bond. */
	void CorrectOverlap_3(const RaggedArray2DT<int>& neighbors, const dArray2DT& coords, double smoothing, double k2, int nip);

	/** solve bond densities one shell at a time */
	void CorrectOverlap_4(const RaggedArray2DT<int>& neighbors, const dArray2DT& coords, double smoothing, double k2, double k_r, int nip);

	/** enforce zero bond density in projected cells */
	void DeactivateFollowerCells(void);
	/*@}*/

	/** (re-)set the equation number for the given group */
	virtual void SetEquationSystem(int group, int start_eq_shift = 0);

	/** set the reference error for the given group */
	void SetReferenceError(int group, double error) const;

	/** return the internal forces for the given solver group associated with the
	 * most recent call to FEManagerT_bridging::FormRHS. */
	const dArray2DT& InternalForce(int group) const;

	/** return the properties map for the given element group. The element group must be
	 * a particle type; otherwise, an exception will be thrown. */
	nMatrixT<int>& PropertiesMap(int element_group);

	/** calculate EAM total electron density at ghost atoms */
	void ElecDensity(int length, dArray2DT& elecdens, dArray2DT& embforce);

	/** add external electron density contribution to ghost atoms */
	void SetExternalElecDensity(const dArray2DT& elecdens, const iArrayT& ghostatoms);
	
	/** add external embedding force contribution to ghost atoms */
	void SetExternalEmbedForce(const dArray2DT& embforce, const iArrayT& ghostatoms);

	/** return the given instance of the ParticlePairT element group or NULL if not found */
	const ParticlePairT* ParticlePair(int instance = 0) const;

	/** \name writing output */
	/*@{*/
	/** initiate the process of writing output from all output sets 
	 * \param time time label associated with the output data */
	virtual void WriteOutput(double time);

	/** write results for a single output set */
	virtual void WriteOutput(int ID, const dArray2DT& n_values, const dArray2DT& e_values) const;

	/* write a snapshot */
	virtual void WriteOutput(const StringT& file, const dArray2DT& coords, const iArrayT& node_map,
		const dArray2DT& values, const ArrayT<StringT>& labels) const;
	/*@}*/

	/** the bridging scale element group */
	BridgingScaleT& BridgingScale(void) const;

protected:

	/** initialize solver information */
	virtual void SetSolver(void);

	/** map coordinates into elements. Temporarily limited to elements
	 * within a single element block */
	void MaptoCells(const iArrayT& nodes, const dArray2DT& coords, iArrayT& cell_num,
		dArray2DT& cell_coords) const;

	/** the EAMT element group */
	EAMT& EAM(void) const;

	/** \name method needed to correct atomistic/continuum overlap */
	/*@{*/
	/** collect nodes and cells in the overlap region. This method collects only free
	 * nodes with ghost atoms in their supports. */
	void CollectOverlapRegion_free(iArrayT& overlap_cell, iArrayT& overlap_node) const;

	/** collect nodes and cells in the overlap region. This method collects all nodes
	 * with ghost nodes in their supports, including both free and projected nodes. */
	//void CollectOverlapRegion_all(iArrayT& overlap_cell, iArrayT& overlap_node) const;

	/** collect bonds between real points and follower points */
	void GhostNodeBonds(const RaggedArray2DT<int>& neighbors, RaggedArray2DT<int>& ghost_neighbors, 
		InverseMapT& overlap_cell_map) const;

	void GhostNodeBonds(const RaggedArray2DT<int>& neighbors,
		RaggedArray2DT<int>& ghost_neighbors, iArrayT& overlap_cell) const;
	
	void GhostNodeBonds_2(const dArrayT& R_i, const dArray2DT& point_coords, 
		const RaggedArray2DT<int>& ghost_neighbors_all, RaggedArray2DT<int>& ghost_neighbors_i, 
		AutoArrayT<int>& overlap_cell_i, AutoArrayT<int>& overlap_node_i) const;

	void GhostNodeBonds_4(double R_i, const dArray2DT& point_coords, 
	const RaggedArray2DT<int>& ghost_neighbors_all, RaggedArray2DT<int>& ghost_neighbors_i, 
	AutoArrayT<int>& overlap_cell_i, AutoArrayT<int>& overlap_node_i) const;

	/** collect list of cells not containing any active bonds */
	void BondFreeElements(const RaggedArray2DT<int>& ghost_neighbors_i, AutoArrayT<int>& bondfree_cell_i) const;

	/** count number of bonds terminating with the domain of each element integration point */
	void CountIPBonds(const ArrayT<int>& elements, iArray2DT& ip_counts) const;

	/** generate "inverse" connectivities for active elements in the support of active nodes.
	 * For each node in the forward map, collect the elements which are part of its support
	 * \param element_group source of element connectivities
	 * \param active_nodes list of nodes for which to generate an element list
	 * \param active_elements list of elements to include in the connectivities
	 * \param transpose_connects returns with the list of active elements in the support of 
	 *        active nodes. The major dimension is the number of active nodes */
	void TransposeConnects(const ContinuumElementT& element_group, const ArrayT<int>& active_nodes,
		const ArrayT<int>& active_elements, RaggedArray2DT<int>& transpose_connects);

	/** compute bond contribution to the nodal internal force
	 * \param R_i bond vector
	 * \param ghost_neighbors neighbor list containing only bond from free points to
	 *        follower points. The follower points only appear in the neighbor lists of the
	 *        free points. The followers do not have neighbors in this list.
	 * \param point_coords global coordinates of all points
	 * \param overlap_node_map map of global node number to index in sum_R_N for nodes in the
	 *        overlap region
	 * \param sum_R_N returns with the bond contribution to all the nodes in the overlap region
	 * \param overlap_cell_i list of cells containing bonds to ghost nodes
	 * \param overlap_atom_i list of atoms with bonds to ghost nodes
	 */
	void ComputeSum_signR_Na(const dArrayT& R_i, const RaggedArray2DT<int>& ghost_neighbors, 
		const dArray2DT& point_coords, const InverseMapT& overlap_node_map, dArrayT& sum_R_N,
		AutoArrayT<int>& overlap_cell_i) const;

	void ComputeSum_signR_Na(const dArrayT& R_i, const RaggedArray2DT<int>& ghost_neighbors, 
		const dArray2DT& point_coords, const InverseMapT& overlap_node_map, dArrayT& sum_R_N) const;

	void ComputeSum_signR_Na_4(double R_i, const RaggedArray2DT<int>& ghost_neighbors, 
		const dArray2DT& point_coords, const InverseMapT& overlap_node_map, dArray2DT& sum_R_N) const;

	/** compute Cauchy-Born contribution to the nodal internal force
	 * \param R bond vector
	 * \param V_0 Cauchy-Born reference volume
	 * \param coarse element group defining shape functions in the overlap region
	 * \param overlap_cell list of elements from coarse in the overlap region
	 * \param overlap_node_map map of global node number to index in sum_R_N for nodes in the
	 *        overlap region
	 */
	void Compute_df_dp_1(const dArrayT& R, double V_0, const ContinuumElementT& coarse, 
		const ArrayT<int>& overlap_cell, const InverseMapT& overlap_node_map, const dArray2DT& rho, 
		dArrayT& f_a, double smoothing, double k2, dArray2DT& df_dp, LAdMatrixT& ddf_dpdp) const;

	void Compute_df_dp(const dArrayT& R, double V_0,
		const ArrayT<char>& cell_type, const InverseMapT& overlap_cell_map, const InverseMapT& overlap_node_map,
		const dArray2DT& rho, dArrayT& f_a, double smoothing, double k2, dArray2DT& df_dp, 
		LAdMatrixT& ddf_dpdp) const;

	void Compute_df_dp_2(const dArrayT& R, double V_0, const ArrayT<char>& cell_type, 
		const InverseMapT& overlap_cell_map, const ArrayT<int>& overlap_node, const InverseMapT& overlap_node_map,
		const iArray2DT& cell_eq_active_i,
		const RaggedArray2DT<int>& inv_connects_i, const RaggedArray2DT<int>& inv_equations_all_i, const RaggedArray2DT<int>& inv_equations_active_i,
		const dArray2DT& rho, dArrayT& f_a, double smoothing, double k2, dArray2DT& df_dp, 
		GlobalMatrixT& ddf_dpdp) const;

	void Compute_df_dp_4(const dArray2DT& shell_bonds, double V_0, const ArrayT<char>& cell_type, 
		const InverseMapT& overlap_cell_map, const ArrayT<int>& overlap_node, const InverseMapT& overlap_node_map,
		const iArray2DT& cell_eq_active_i,
		const RaggedArray2DT<int>& inv_connects_i, const RaggedArray2DT<int>& inv_equations_i,
		const dArray2DT& rho, dArray2DT& f_a, double smoothing, double k2, dArray2DT& df_dp, 
		GlobalMatrixT& ddf_dpdp) const;
	/*@}*/

	/** group bonds into shells */
	int NumberShells(const dArray2DT& bonds, iArrayT& shell, dArrayT& shell_bond_length) const;

private:

	/** \name ghost node information */
	/*@{*/
	/** list of my ghost nodes */
	iArrayT fGhostNodes;

	/** list of my non-ghost nodes */
	iArrayT fNonGhostNodes;
	
	/** ghost nodes pseudo-equations */
	iArray2DT fGhostNodesEquations;
	
	/** map from ghost node id to row in FEManagerT_bridging::fGhostNodesEquations */
	InverseMapT fGhostIdToIndex;
	/*@}*/

	/** projection/interpolation operator */
	BridgingScaleT* fBridgingScale;
	
	/** EAMFCC class */
	EAMFCC3D* fEAMFCC3D;
	
	/** EAMT class */
	EAMT* fEAMT;
	
	/** \name follower node information */
	/*@{*/
	/** map data of follower points into the mesh */
	PointInCellDataT fFollowerCellData;
	/*@}*/
	
	/** \name the driven solution */
	/*@{*/
	/** the KBC_ControllerT applying the external solution */
	KBC_ControllerT* fSolutionDriver;
	
	/** map data of driver points into the mesh */
	PointInCellDataT fDrivenCellData;
	
	/** list of projected nodes */
	iArrayT fProjectedNodes;
	
	/** projected solution */
	dArray2DT fProjection;
	
	/** cumulative update for each solver group */
	ArrayT<dArrayT> fCumulativeUpdate;

	/** active element map. Skip calculation of element without active nodes
	 * and not participating in the overlap layer. */
	ArrayT<ElementCardT::StatusT> fElementStatus;
	/*@}*/
	
	/** \name external force vector by group */
	/*@{*/
	ArrayT<const dArrayT*> fExternalForce;	
	ArrayT<const dArray2DT*> fExternalForce2D;
	ArrayT<const iArrayT*>   fExternalForce2DNodes;
	ArrayT<iArray2DT>        fExternalForce2DEquations;
	/*@}*/
};

/* inlines */

/* set pointer to an external force vector or pass NULL to clear */
inline void FEManagerT_bridging::SetExternalForce(int group, const dArrayT& external_force)
{
	fExternalForce[group] = &external_force;
}

// clear pointers to external (coupling) force
inline void FEManagerT_bridging::ClearExternalForce(int group)
{
	fExternalForce[group] = NULL;
	fExternalForce2D[group] = NULL;
}

} /* namespace Tahoe */

#endif /* BRIDGING_ELEMENT */
#endif /* _FE_MANAGER_BRIDGING_H_ */
