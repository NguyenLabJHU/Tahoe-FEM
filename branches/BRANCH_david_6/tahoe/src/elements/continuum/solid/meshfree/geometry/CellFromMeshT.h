/* $Id: CellFromMeshT.h,v 1.4 2005-09-29 19:14:47 jcmach Exp $ */
#ifndef _CELL_FROM_MESH_T_H_
#define _CELL_FROM_MESH_T_H_

/* base class */
#include "CellGeometryT.h"

/* direct members */
#include "iArray2DT.h"

namespace Tahoe {

/** class which computes the nodal boundary integrals using surfaces of predefined
 * sub-volumes of a mesh used to define the body geometry. */
class CellFromMeshT: public CellGeometryT
{
public:

	/** constructor */
	CellFromMeshT(const ElementSupportT& support, bool isAxisymmetric);
	CellFromMeshT(void);

	/** echo element connectivity data. Reads parameters that define
	 * which nodes belong to this ParticleT group. */
	virtual void DefineElements(const ArrayT<StringT>& block_ID, const ArrayT<int>& mat_index);
	
	/** generate data structures for integration over the body boundary */
	virtual void BoundaryShapeFunctions(RaggedArray2DT<double>& phis, RaggedArray2DT<int>& supports, dArray2DT& normals);	
	
	/** compute B matrices for strain smoothing/nodal integration */
	virtual void ComputeBMatrices(RaggedArray2DT<int>& nodalCellSupports, RaggedArray2DT<dArrayT>& bVectorArray,
				      dArrayT& cellVolumes, dArray2DT& cellCentroids, RaggedArray2DT<double>& circumferential_B);

	/** compute Bprime matrices for strain smoothing/nodal integration */
 	virtual void ComputeBprimeMatricesSS(RaggedArray2DT<dMatrixT>& bprimeVectorArray, const RaggedArray2DT<int>& cellSupports,
					   const RaggedArray2DT<dArrayT>& bVectorArray, const dArrayT& cellVolumes, const dArray2DT& cellCentroids,
					   dArray2DT& Ymatrices); 

	
protected:

	/** block ID's defining the meshless domain */
	ArrayT<StringT> fBlockID;
};

} /* namespace Tahoe */

#endif /* _CELL_FROM_MESH_T_H_ */


