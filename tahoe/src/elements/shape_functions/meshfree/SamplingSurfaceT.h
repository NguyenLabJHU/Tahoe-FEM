/* $Id: SamplingSurfaceT.h,v 1.3 2002-07-05 22:28:37 paklein Exp $ */
/* created: paklein (10/19/2000)                                          */

#ifndef _SAMPLING_SURFACE_T_H_
#define _SAMPLING_SURFACE_T_H_

/* direct members */
#include "GeometryT.h"
#include "iArrayT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"
#include "RaggedArray2DT.h"
#include "LocalArrayT.h"
#include "dMatrixT.h"

namespace Tahoe {

/* forward declarations */
class MeshFreeSupportT;
template <class TYPE> class ArrayT;
class SurfaceShapeT;

class SamplingSurfaceT
{
public:

	/* constructor */
	SamplingSurfaceT(GeometryT::CodeT code, int num_facet_nodes,
		int num_samples, MeshFreeSupportT& mf_support);

	/* destructor */
	~SamplingSurfaceT(void);

	/* dimensions */
	int NumberOfFacets(void) const;
	int SamplesPerFacet(void) const;
	int NodesPerFacet(void) const;
	
	/* compute/store shape functions at sampling points on surface */
	void SetSamplingPoints(const dArray2DT& facet_coords,
		const iArray2DT& facet_connects);

	/* reset database for facets affected by given nodes */
	void ResetSamplingPoints(const iArrayT& changed_nodes);

	/* retrieve stored values */
	void LoadSamplingPoint(int facet, int point, iArrayT& nodes, dArrayT& phi,
		dArray2DT& Dphi);

	/* reference surface configuration, where Q brings a vector to the
	 * local frame by:
	 *
	 *       t_local = Q^T t_global
	 *
	 * normal is returned in the global frame.
	 */
	void ReferenceConfiguration(int facet, int point, dMatrixT& Q, dArrayT& normal);
	
	/* facet coordinates */
	void FacetCoordinates(int facet, dArray2DT& coordinates) const;

	/* user-definable flag per facet */
	int ChangeFlags(int from, int to); // return the number changed
	void SetFlags(int value);
	int& Flag(int facet);

	/* write storage statistics */
	void WriteStatistics(ostream& out) const;

private:

	/* compute/store shape function data - for set neighbors,
	 * (facets != NULL) => selective recomputation */
	void SetFieldData(const ArrayT<int>* facets);
	
protected:
	
	/* parameters */
	GeometryT::CodeT fCode;
	int fNumFacetNodes;
	int fNumSamples;

	/* geometry data */
	dArray2DT fFacetCoords;
	iArray2DT fFacetConnects;

	/* meshfree support */
	MeshFreeSupportT& fMFSupport;

	/* surface geometry */
	LocalArrayT fLocFacetCoords;
	SurfaceShapeT* fSurfaceShape;
		
	/* shape function database */
	int fNumFacets;
	RaggedArray2DT<int> fNeighborData; // [nfacet x nip] x [nnd]
	RaggedArray2DT<double> fPhiData;   // [nfacet x nip] x [nnd]
	RaggedArray2DT<double> fDPhiData;  // [nfacet x nip] x [nsd x nnd]
	
	/* user-definable flag per facet */
	iArrayT fFlag;

	/* work space */
	dMatrixT fQ; // need for extracting surface normal
	
	/* sampling points to skip -> parallel execution */
	iArrayT fSkipPoints;
	// not implemented
};

/* inlines */
inline int SamplingSurfaceT::ChangeFlags(int from, int to) { return fFlag.ChangeValue(from, to); }
inline void SamplingSurfaceT::SetFlags(int value) { fFlag = value; }
inline int& SamplingSurfaceT::Flag(int facet) { return fFlag[facet]; }

/* dimensions */
inline int SamplingSurfaceT::NumberOfFacets(void) const { return fNumFacets; }
inline int SamplingSurfaceT::SamplesPerFacet(void) const { return fNumSamples; }
inline int SamplingSurfaceT::NodesPerFacet(void) const { return fNumFacetNodes; }

} // namespace Tahoe 
#endif /* _SAMPLING_SURFACE_T_H_ */
