/* $Id: FEManagerT_bridging.h,v 1.11.12.2 2004-06-07 13:50:51 paklein Exp $ */
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

namespace Tahoe {

/* forward declarations */
class ParticleT;
class BridgingScaleT;
class KBC_ControllerT;
class dSPMatrixT;

/** extension of FEManagerT for bridging scale calculations */
class FEManagerT_bridging: public FEManagerT
{
public:

	/** constructor */
	FEManagerT_bridging(ifstreamT& input, ofstreamT& output, CommunicatorT& comm,
		ifstreamT& bridging_input);

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

	/** set pointer to an external force vector for the given field */
	void SetExternalForce(const StringT& field, const dArray2DT& external_force, const iArrayT& activefenodes);
	/*@}*/

	/** \name ghost nodes 
	 * The ghost node database must be initialized by calling
	 * FEManagerT_bridging::InitGhostNodes before accessing the lists.*/
	/*@{*/
	/** initialize the ghost node information 
	 * \param include_image_nodes flag to indicate whether image nodes should be
	 *        included in the list of non-ghost nodes */
	void InitGhostNodes(bool include_image_nodes);

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
	void InitInterpolation(const iArrayT& nodes, const StringT& field,
		NodeManagerT& node_manager);

	/** field interpolations. Interpolate the field to the nodes initialized
	 * with the latest call to FEManagerT_bridging::InitInterpolation. */
	void InterpolateField(const StringT& field, int order, dArray2DT& nodal_values);

	/** return the interpolation matrix associated with the active degrees
	 * of freedom */
	void InterpolationMatrix(const StringT& field, dSPMatrixT& G_Interpolation) const;

	/** compute global interpolation matrix for all nodes whose support intersects the MD 
	 *  region, i.e. N_{I}(X_{\alpha}) */
	void Ntf(dSPMatrixT& ntf, const iArrayT& atoms, iArrayT& activefenodes) const;

	/** initialize projection data. Initialize data structures needed to project
	 * field values to the given list of points. Requires that this FEManagerT has
	 * a BridgingScaleT in its element list. */
	void InitProjection(CommManagerT& comm, const iArrayT& nodes, const StringT& field,
		NodeManagerT& node_manager, bool make_inactive);

	/** indicate whether image nodes should be included in the projection */
	virtual bool ProjectImagePoints(void) const;

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
	/*@}*/

	/** (re-)set the equation number for the given group */
	virtual void SetEquationSystem(int group);

	/** set the reference error for the given group */
	void SetReferenceError(int group, double error) const;

	/** return the internal forces for the given solver group associated with the
	 * most recent call to FEManagerT_bridging::FormRHS. */
	const dArray2DT& InternalForce(int group) const;

	/** return the properties map for the given element group. The element group must be
	 * a particle type; otherwise, an exception will be thrown. */
	nMatrixT<int>& PropertiesMap(int element_group);

protected:

	/** initialize solver information */
	virtual void SetSolver(void);

	/** map coordinates into elements. Temporarily limited to elements
	 * within a single element block */
	void MaptoCells(const iArrayT& nodes, const dArray2DT& coords, iArrayT& cell_num,
		dArray2DT& cell_coords) const;

	/** the bridging scale element group */
	BridgingScaleT& BridgingScale(void) const;

private:

	/** input parameters for bridging parameters */
	ifstreamT& fBridgingIn;

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
	
	/** projected solution */
	dArray2DT fProjection;
	
	/** cumulative update for each solver group */
	ArrayT<dArrayT> fCumulativeUpdate;
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

} /* namespace Tahoe */

#endif /* BRIDGING_ELEMENT */
#endif /* _FE_MANAGER_BRIDGING_H_ */
