/* $Id: CellFromMeshT.h,v 1.2 2005-01-25 18:22:11 cjkimme Exp $ */
#ifndef _CELL_FROM_MESH_T_H_
#define _CELL_FROM_MESH_T_H_

/* base class */
#include "CellGeometryT.h"

/* direct members */
#include "iArray2DT.h"

namespace Tahoe {

/** base class for particle types */
class CellFromMeshT: public CellGeometryT
{
public:

	/** constructor */
	CellFromMeshT(const ElementSupportT& support, bool isAxisymmetric);
	CellFromMeshT(void);

	/** destructor */
	~CellFromMeshT(void);
	
	/** echo element connectivity data. Reads parameters that define
	 * which nodes belong to this ParticleT group. */
	virtual void DefineElements(const ArrayT<StringT>& block_ID, const ArrayT<int>& mat_index);
	
	/** generate data structures for integration over the body boundary */
	virtual void BoundaryShapeFunctions(RaggedArray2DT<double>& phis, RaggedArray2DT<int>& supports, dArray2DT& normals);	
	
	/** compute B matrices for strain smoothing/nodal integration */
	virtual void ComputeBMatrices(RaggedArray2DT<int>& nodalCellSupports, RaggedArray2DT<dArrayT>& bVectorArray,
									dArrayT& cellVolumes, dArray2DT& cellCentroids, RaggedArray2DT<double>& circumferential_B);
	
protected:
	
	/** connectivities used to compute B matrices. */
	iArray2DT fElementConnectivities;

};

} /* namespace Tahoe */

#endif /* _CELL_FROM_MESH_T_H_ */


