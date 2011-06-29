/* $Id: MeshFreeSurfaceSupportT.h,v 1.4 2002-10-20 22:49:41 paklein Exp $ */
/* created: paklein (02/22/2000)                                          */

#ifndef _MF_SURFACE_SUPPORT_T_H_
#define _MF_SURFACE_SUPPORT_T_H_

#include "Environment.h"

/* direct members */
#include "iArrayT.h"
#include "dArray2DT.h"
#include "nVariArray2DT.h"
#include "VariRaggedArray2DT.h"
#include "iAutoArrayT.h"

namespace Tahoe {

/* forward declarations */
class MeshFreeSupportT;
class dArray2DT;
class SurfaceShapeT;
class LocalArrayT;

class MeshFreeSurfaceSupportT
{
public:

	/* constructor */
	MeshFreeSurfaceSupportT(MeshFreeSupportT& mf_support,
		SurfaceShapeT& ref_surface_shape, LocalArrayT& ref_loc_coords,
		const dArray2DT& facet_coords, int num_facet_nodes,
		bool storeshape);

	/* (re-)compute shape functions - NULL resets all facets */
	const dArray2DT& FacetCoords(void) const;
	void ResetFacets(const ArrayT<int>* reset_facets = NULL);

	/* "load" data for the specified facet for all integration points */
//TEMP
	void LoadData(int facet,
		iArrayT& neighbors_1, dArray2DT& phi_1, ArrayT<dArray2DT>& Dphi_1,
		iArrayT& neighbors_2, dArray2DT& phi_2, ArrayT<dArray2DT>& Dphi_2);

	int NumberOfNeighbors(int facet) const;
	void LoadData(int facet, int side, iArrayT& neighbors, dArray2DT& phi,
		ArrayT<dArray2DT>& Dphi);

	/* access to neighbors database */
	const ArrayT<int>& NeighborCounts(int side) const;
	const RaggedArray2DT<int>& Neighbors(int side) const;

	/* list of nodes used in the connectivities */
//	const iArrayT& NodesUsed(void) const;

	/* write MLS statistics */
	void WriteStatistics(ostream& out) const;

private:

	/* make inverse map: inv_map[map[i] - shift] = i */
	void MakeInverseMap(const iAutoArrayT& map, AutoArrayT<int>& inv_map,
		int& shift) const;

	/* stores data for the given facet/side */
	void StoreData(int facet, int side,
		ArrayT< AutoArrayT<int> >& neighbors,
		ArrayT< AutoArrayT<double> >& phi,
		ArrayT< AutoArrayT<double> >& Dphi);

	/* returns the characteristic facet dimensions */
	double FacetSize(const dArray2DT& facet_coords) const;

protected:

	/* parameters */
	bool fStoreShape;

/* meshfree support */
	MeshFreeSupportT& fMFSupport;

	/* surface shape function and their local coordinate array */
	SurfaceShapeT& fRefSurfaceShape;
	LocalArrayT&   fRefLocCoords;

	/* facet coordinate data */
	int fNumFacetNodes;
	const dArray2DT& fFacetCoords;
	
	/* nodal shape function database */
	AutoArrayT<int>            fCount_1;
	VariRaggedArray2DT<int>    fNeighbors_1; // rows: [nip] x [nnd]
	VariRaggedArray2DT<double> fPhi_1;       // rows: [nip] x [nnd]
	VariRaggedArray2DT<double> fDPhi_1;      // rows: [nip] x [nsd] x [nnd]

	AutoArrayT<int>            fCount_2;
	VariRaggedArray2DT<int>    fNeighbors_2;
	VariRaggedArray2DT<double> fPhi_2;
	VariRaggedArray2DT<double> fDPhi_2;
	
private:

	/* database management workspace */
	iAutoArrayT fall_neighbors; // full ip +/- neighborhood
	AutoArrayT<int> fnode_map;  // global->element node map
	AutoArrayT<double> fphi;
	AutoArrayT<double> fDphi;
};

/* inlines */
inline const dArray2DT& MeshFreeSurfaceSupportT::FacetCoords(void) const
{
	return fFacetCoords;
}

/* access to neighbors database */
inline const ArrayT<int>& MeshFreeSurfaceSupportT::NeighborCounts(int side) const
{
#if __option(extended_errorcheck)
	if (side != 0 && side != 1) throw ExceptionT::kOutOfRange;
#endif

	if (side == 0)
		return fCount_1;
	else
		return fCount_2;
}

inline const RaggedArray2DT<int>& MeshFreeSurfaceSupportT::Neighbors(int side) const
{
#if __option(extended_errorcheck)
	if (side != 0 && side != 1) throw ExceptionT::kOutOfRange;
#endif

	if (side == 0)
		return fNeighbors_1;
	else
		return fNeighbors_2;
}

inline int MeshFreeSurfaceSupportT::NumberOfNeighbors(int facet) const
{
	return fCount_1[facet] + fCount_2[facet];
}


/* list of nodes used in the connectivities */
//inline const iArrayT& MeshFreeSurfaceSupportT::NodesUsed(void) const
//{
//	return fNodesUsed;
//}
} // namespace Tahoe 
#endif /* _MF_SURFACE_SUPPORT_T_H_ */
