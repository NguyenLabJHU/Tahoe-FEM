/* $Id: MeshfreeBridgingT.h,v 1.2.12.1 2003-10-16 12:51:29 paklein Exp $ */
#ifndef _MESHFREE_BRIDGING_SCALE_T_H_
#define _MESHFREE_BRIDGING_SCALE_T_H_

/* base class */
#include "BridgingScaleT.h"

namespace Tahoe {

/* forward declarations */
class MLSSolverT;

/** bridging class that uses moving least squares fits to transfer information
* from arbitray points onto the nodes. */
class MeshfreeBridgingT: public BridgingScaleT
{
public:

	/** constructor */
	MeshfreeBridgingT(const ElementSupportT& support, const FieldT& field,
		const SolidElementT& solid);

	/** destructor */
	~MeshfreeBridgingT(void);

	/** \name "project" external field onto mesh */
	/*@{*/
	/** initialize projection data. Maps points into the initial or current 
	 * configuration of the mesh. Only one of the two pointers must be passed as NULL.
	 * \param points_used indecies of points in the coordinate lists which should be
	 *        mapped into cells
	 * \param init_coords point to the initial coordinates array if available
	 * \param curr_coords point to the initial coordinates array if available
	 * \param point_in_cell destination for map data
	 */
	virtual void InitProjection(const iArrayT& points_used, const dArray2DT* init_coords, 
		const dArray2DT* curr_coords, PointInCellDataT& cell_data);

	/** project the point values onto the mesh. Requires a previous call to
	 * BridgingScaleT::InitProjection to initialize the PointInCellDataT */
	virtual void ProjectField(const PointInCellDataT& cell_data,
		const dArray2DT& point_values, dArray2DT& projection);
	/*@}*/

protected:

	/** determines points in the neighborhoods of nodes of each non-empty cell.
	 * Only one of the two pointers must be passed as NULL. This method also calls
	 * MeshfreeBridgingT::MaptoCells to determine the non-empty cells and the list
	 * of nodes in those cells.
	 * \param points_used indecies of points in the coordinate lists which should be
	 *        mapped into cells
	 * \param init_coords point to the initial coordinates array if available
	 * \param curr_coords point to the initial coordinates array if available
	 * \param point_in_cell destination for map data
	 */
	void BuildNodalNeighborhoods(const iArrayT& points_used, const dArray2DT* init_coords, 
		const dArray2DT* curr_coords, PointInCellDataT& cell_data);

private:

	/** moving least squares solver */
	MLSSolverT* fMLS;
	
	/** support size of each atom in the crystal, not just the points in cells */
	dArray2DT fSupport;
};

} /* namespace Tahoe */

#endif /* _MESHFREE_BRIDGING_SCALE_T_H_ */
