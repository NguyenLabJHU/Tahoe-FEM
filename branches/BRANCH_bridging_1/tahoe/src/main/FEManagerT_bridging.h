/* $Id: FEManagerT_bridging.h,v 1.1.2.8 2003-02-23 02:41:23 paklein Exp $ */
#ifndef _FE_MANAGER_BRIDGING_H_
#define _FE_MANAGER_BRIDGING_H_

/* base class */
#include "FEManagerT.h"

/* direct members */
#include "PointInCellDataT.h"

namespace Tahoe {

/* forward declarations */
class ParticleT;
class BridgingScaleT;
class KBC_PrescribedT;
class dSPMatrixT;

class FEManagerT_bridging: public FEManagerT
{
public:

	/** constructor */
	FEManagerT_bridging(ifstreamT& input, ofstreamT& output, CommunicatorT& comm,
		ifstreamT& bridging_input);

	/** \name ghost nodes 
	 * The ghost node database must be initialized by calling
	 * FEManagerT_bridging::InitGhostNodes before accessing the lists.*/
	/*@{*/
	/** initialize the ghost node information */
	void InitGhostNodes(void);

	/** return list of ghost nodes */
	const iArrayT& GhostNodes(void) const { return fGhostNodes; };

	/** return list of ghost nodes */
	const iArrayT& NonGhostNodes(void) const { return fNonGhostNodes; };

	/** compute the ghost-nonghost part of the stiffness matrix */
	void Form_G_NG_Stiffness(const StringT& field, dSPMatrixT& K_G_NG);
	/*@}*/

	/** write field values for the given nodes */
	void SetFieldValues(const StringT& field, const iArrayT& nodes, 
		const dArray2DT& values);

	/** \name interpolation and projection operators */
	/*@{*/
	/** initialize interpolation data. Initialize data structures needed to interpolate
	 * field values to the given list of points. Requires that this FEManagerT has
	 * a BridgingScaleT in its element list. */
	void InitInterpolation(const iArrayT& nodes, const StringT& field,
		NodeManagerT& node_manager);

	/** field interpolations. Interpolate the field to the nodes initialized
	 * with the latest call to FEManagerT_bridging::InitInterpolation. */
	void InterpolateField(const StringT& field, dArray2DT& nodal_values);

	/** return the interpolation matrix associated with the active degrees
	 * of freedom */
	void InterpolationMatrix(const StringT& field, dSPMatrixT& G_Interpolation) const;

	/** initialize projection data. Initialize data structures needed to project
	 * field values to the given list of points. Requires that this FEManagerT has
	 * a BridgingScaleT in its element list. */
	void InitProjection(const iArrayT& nodes, const StringT& field,
		NodeManagerT& node_manager);

	/** project the point values onto the mesh. Project to the nodes using
	 * projection initialized with the latest call to FEManagerT_bridging::InitProjection. */
	void ProjectField(const StringT& field, NodeManagerT& node_manager);
	/*@}*/

	/** \name solver control */
	/*@{*/
	/** the residual for the given group. The array contains the residual from
	 * the latest call to FEManagerT::FormRHS */
	const dArrayT& Residual(int group) const;

	/** set the reference error for the given group */
	void SetReferenceError(int group, double error) const;
	/*@}*/

private:

	/** map coordinates into elements. Temporarily limited to elements
	 * within a single element block */
	void MaptoCells(const iArrayT& nodes, const dArray2DT& coords, iArrayT& cell_num,
		dArray2DT& cell_coords) const;

	/** the particle element group */
	ParticleT& Particle(void) const;

	/** the bridging scale element group */
	BridgingScaleT& BridgingScale(void) const;

private:

	/** input parameters for bridging parameters */
	ifstreamT& fBridgingIn;

	/** \name ghost node information */
	/*@{*/
	ParticleT* fParticle;
	
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
	KBC_PrescribedT* fSolutionDriver;
	
	/** map data of driver points into the mesh */
	PointInCellDataT fDrivenCellData;
	
	/** projected solution */
	dArray2DT fProjection;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _FE_MANAGER_BRIDGING_H_ */
