/* $Id: FEManagerT_bridging.h,v 1.1.2.2 2003-02-10 02:22:04 paklein Exp $ */
#ifndef _FE_MANAGER_BRIDGING_H_
#define _FE_MANAGER_BRIDGING_H_

/* base class */
#include "FEManagerT.h"

/* direct members */
#include "PointInCellDataT.h"

namespace Tahoe {

/* forward declarations */
class BridgingScaleT;
class KBC_PrescribedT;

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
	/*@}*/

	/** initialize interpolation data. Initialize data structures needed to interpolate
	 * field values to the given list of points. Requires that this FEManagerT has
	 * a BridgingScaleT in its element list. */
	void InitInterpolation(const iArrayT& nodes, const StringT& field,
		NodeManagerT& node_manager);

	/** initialize projection data. Initialize data structures needed to project
	 * field values to the given list of points. Requires that this FEManagerT has
	 * a BridgingScaleT in its element list. */
	void InitProjection(const iArrayT& nodes, const StringT& field,
		NodeManagerT& node_manager);

private:

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
	/*@}*/
	
	/** \name follower node information */
	/*@{*/
	BridgingScaleT* fBridgingScale;
	
	/** map data of follower points into the mesh */
	PointInCellDataT fFollowerCellData;
	/*@}*/
	
	/** \name the driven solution */
	/*@{*/
	/** the KBC_ControllerT applying the external solution */
	KBC_PrescribedT* fSolutionDriver;
	
	/** list of nodes in elements being driven */
	iArrayT fDrivenCellNodes;
	
	/** map data of driver points into the mesh */
	PointInCellDataT fDrivenCellData;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _FE_MANAGER_BRIDGING_H_ */
