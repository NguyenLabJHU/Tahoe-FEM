/* $Id: MeshfreeBridgingT.h,v 1.9 2005-04-16 01:59:37 paklein Exp $ */
#ifndef _MESHFREE_BRIDGING_SCALE_T_H_
#define _MESHFREE_BRIDGING_SCALE_T_H_

/* base class */
#include "BridgingScaleT.h"

/* direct members */
#include "MeshFreeSupport3DT.h"

namespace Tahoe {

/* forward declarations */
class MLSSolverT;

/** bridging class that uses moving least squares fits to transfer information
* from arbitray points onto the nodes. */
class MeshfreeBridgingT: public BridgingScaleT
{
public:

	/** constructor */
	MeshfreeBridgingT(const ElementSupportT& support);

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
	virtual void InitProjection(CommManagerT& comm, const iArrayT& points_used, const dArray2DT* init_coords, 
		const dArray2DT* curr_coords, PointInCellDataT& cell_data);

	/** project the point values onto the mesh. Requires a previous call to
	 * BridgingScaleT::InitProjection to initialize the PointInCellDataT */
	virtual void ProjectField(const PointInCellDataT& cell_data,
		const dArray2DT& point_values, dArray2DT& projection);

	/** compute the coarse scale part of the source field */
	virtual void CoarseField(const PointInCellDataT& cell_data,
		const dArray2DT& field, dArray2DT& coarse) const;

	/** indicate whether image nodes should be included in the projection */
	virtual bool ProjectImagePoints(void) const { return false; };

	/** collect the cells without any free nodes */
	virtual void CollectProjectedCells(const PointInCellDataT& cell_data, iArrayT& cells) const;

	/** return list of projected nodes */
	virtual void CollectProjectedNodes(const PointInCellDataT& cell_data, iArrayT& nodes) const;

	/** compute \f$ B_{\hat{U}U} \f$. Compute the matrix defining the projection of data from free cell
	 * nodes to driven cell nodes which arises because some projection source points lie in cells containing
	 * both free and projection nodes
	 * \param projection projection data containing both information about the projection from source points
	 *        to driven cell nodes \f$ B_{\hat{U}Q} \f$ and interpolation from free cell nodes to the projection 
	 *        source points \f$ N_{QU} \f$.
	 * \param B_hatU_U destination for the projection matrix \f$ B_{\hat{U}U} \f$ */
	virtual void Compute_B_hatU_U(const PointInCellDataT& projection, InterpolationDataT& B_hatU_U) const;
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;
	
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** determines points in the neighborhoods of nodes of each non-empty cell.
	 * Only one of the two pointers must be passed as NULL. This method also calls
	 * MeshfreeBridgingT::MaptoCells to determine the non-empty cells and the list
	 * of nodes in those cells.
	 * \param points_used indecies of points in the coordinate lists which should be
	 *        mapped into cells
	 * \param init_coords pointer to the initial coordinates array if available
	 * \param curr_coords pointer to the initial coordinates array if available
	 * \param point_in_cell destination for map data
	 */
	void BuildNodalNeighborhoods(CommManagerT& comm, const iArrayT& points_used, const dArray2DT* init_coords, 
		const dArray2DT* curr_coords, PointInCellDataT& cell_data);

	/** determines points in the neighborhoods of points used. This method must be
	 * called after MeshfreeBridgingT::fSupport has been found, which occurs during the
	 * call to MeshfreeBridgingT::BuildNodalNeighborhoods.
	 * \param points_used indecies of points in the coordinate lists which should be
	 *        mapped into cells
	 * \param point_coords point coordinates
	 * \param point_in_cell destination for map data
	 */
	void BuildPointNeighborhoods(const iArrayT& points_used, const dArray2DT& point_coords, 
		PointInCellDataT& cell_data);

private:

	/** moving least squares solver */
	MLSSolverT* fMLS;
	
	/** support size parameters of each atom in the crystal, not just the points in cells */
	dArray2DT fSupportParams;

	/** needed to define MLSSolverT parameters */
	MeshFreeSupport3DT fMeshFreeSupport;
};

} /* namespace Tahoe */

#endif /* _MESHFREE_BRIDGING_SCALE_T_H_ */
