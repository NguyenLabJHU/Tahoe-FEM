/* $Id: InterpolationDataT.h,v 1.3 2004-06-26 05:53:19 paklein Exp $ */
#ifndef _INTERPOLATION_DATA_T_H_
#define _INTERPOLATION_DATA_T_H_

/* direct members */
#include "RaggedArray2DT.h"
#include "InverseMapT.h"

namespace Tahoe {

/* forward declarations */
class iArray2DT;
class dArray2DT;
class dArrayT;

/** collection of information for interpolating data */
class InterpolationDataT
{
public:

	/** constructor */
	InterpolationDataT(void) { };

	/** \name accessors */
	/*@{*/
	/** neighbors for each entry in InterpolationDataT::Map */
	RaggedArray2DT<int>& Neighbors(void) { return fNeighbors; };
	const RaggedArray2DT<int>& Neighbors(void) const { return fNeighbors; };

	/** interpolation weights */
	RaggedArray2DT<double>& NeighborWeights(void) { return fNeighborWeights; };
	const RaggedArray2DT<double>& NeighborWeights(void) const { return fNeighborWeights; };

	/** map of global number to corresponding row in InterpolationDataT::Neighbors
	 * and InterpolationDataT::NeighborWeights */
	InverseMapT& Map(void) { return fMap; };
	const InverseMapT& Map(void) const { return fMap; };
	/*@}*/

	/** \name transpose the given interpolation data 
	 * \param map map from global id of interpolation point to the row in the
	 *        interpolation data.
	 * \param neighbors ids of the neighbors for each interpolation point
	 * \param neighbor_weights interpolation weights of the neighbors for each 
	 *        interpolation point
	 */
	/*@{*/
	/** transpose interpolation from an arbitrary set of nodes */
	void Transpose(const InverseMapT& map, const RaggedArray2DT<int>& neighbors,
		const RaggedArray2DT<double>& neighbor_weights);

	/** transpose interpolation from regular connectivities */
	void Transpose(const InverseMapT& map, const iArray2DT& neighbors,
		const dArray2DT& neighbor_weights);
	/*@}*/

	/** return the interpolation data as a sparse matrix */
	void ToMatrix(iArrayT& r, iArrayT& c, dArrayT& v) const;
	
private:

	/** map of global number to corresponding row in InterpolationDataT::fNeighbors
	 * and InterpolationDataT::fNeighborWeights */
	InverseMapT fMap;
	
	/** neighbors for each entry in InterpolationDataT::fMap */
	RaggedArray2DT<int> fNeighbors;

	/** interpolation weights */
	RaggedArray2DT<double> fNeighborWeights;	
};

} /* namespace Tahoe */

#endif /* _INTERPOLATION_DATA_T_H_ */
