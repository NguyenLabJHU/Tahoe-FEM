/* $Id: Front2DT.cpp,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: paklein (03/18/1999)                                          */

#include "Front2DT.h"

#include "FrontSegmentT.h"
#include "FrontNodeT.h"
#include "iArrayT.h"

/* constants */
const double Pi = acos(-1.0);

/* constructor */
Front2DT::Front2DT(double cone, double da, double da_s, int num_pts):
	FrontT(2, 2, cone, da, da_s, num_pts)
{

}

/* construct initial front */
void Front2DT::Initialize(const dArray2DT& facet_coords, const iArrayT& fr_facets,
		const iArrayT& fr_edges)
{
	/* expecting line segments */
	if (facet_coords.MinorDim() != 2*2) throw eGeneralFail;

	for (int i = 0; i < fr_facets.Length(); i++)
	{
		double* x1 = facet_coords(fr_facets[i]);
		double* x2 = x1 + 2;

		double* x = (fr_edges[i] == 0) ? x1 : x2;
		double v_n[2];
		if (fr_edges[i] == 0)
		{
			x = x1;		
			v_n[0] = x1[0] - x2[0];
			v_n[1] = x1[1] - x2[1];
		}
		else
		{
			x = x2;
			v_n[0] = x2[0] - x1[0];
			v_n[1] = x2[1] - x1[1];
		}

		/* construct node */
		FrontNodeT* node = new FrontNodeT(2, x, v_n, NULL, fcone, fda*fda_s, fnum_pts);
		if (!node) throw eOutOfMemory;
	
		/* store */
		fFrontNodes.Append(node);
	}
}  		

/* extend front at the specified points along the given direction - returns new
* cutting facets */
const dArray2DT& Front2DT::NewFacets(const ArrayT<int>& extend_pts,
	const ArrayT<int>& extend_dir)
{
	/* dimension check */
	if (extend_pts.Length() != extend_dir.Length()) throw eSizeMismatch;

	/* each extension generates a new facets */
	fNewFacetMan.SetMajorDimension(extend_pts.Length(), false);

	/* generate front facets */
	double n_new[2];
	for (int i = 0; i < extend_pts.Length(); i++)
	{
		/* segment coordinates */
		double* x1 = fNewFacets(i);
		double* x2 = x1 + 2;
	
		/* node and sampling point */
		int node = extend_pts[i];
		int dir  = extend_dir[i];
		FrontNodeT& frontnode = *fFrontNodes[node];
	
		/* origin and direction */
		const double* x = frontnode.Coords();
		const double* n = frontnode.Direction(dir);
		
		/* new front facet */
		x1[0] = x[0];
		x1[1] = x[1];
		
		x2[0] = x1[0] + fda*n[0];
		x2[1] = x1[1] + fda*n[1];
	
		/* forward direction */
		n_new[0] = x2[0] - x1[0];
		n_new[1] = x2[1] - x1[1];
	
		/* reset the front node */
		frontnode.Reset(2, x2, n_new, NULL, fcone, fda*fda_s, fnum_pts);
	}

	return fNewFacets;
}
