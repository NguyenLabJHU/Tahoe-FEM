/* $Id: Front3DT.h,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: paklein (03/18/1999)                                          */

#ifndef _FRONT_3D_T_H_
#define _FRONT_3D_T_H_

/* base class */
#include "FrontT.h"

/* forward declarations */
class FrontSegmentT;

class Front3DT: public FrontT
{
public:

	/* constructor */
	Front3DT(double cone, double da, double da_s, int num_pts);

	/* destructor */
	~Front3DT(void);

	/* construct initial front */
	virtual void Initialize(const dArray2DT& facet_coords, const iArrayT& fr_facets,
		const iArrayT& fr_edges);
	
	/* extend front at the specified points along the given direction - returns new
	 * cutting facets */
	virtual const dArray2DT& NewFacets(const ArrayT<int>& extend_pts,
		const ArrayT<int>& extend_dir);

	/* activate kink control by specifying angle > 0 degrees */
	void SetKinkAngle(double angle);

private:

	// OBSOLETE FUNCTION
	const dArray2DT& NewFacets_insert(const ArrayT<int>& extend_pts,
		const ArrayT<int>& extend_dir);

	/* check for minimum included angle - returns facets added to remove
	 * kinks in the front */
	const dArray2DT& RemoveKinks(double min_angle);
		
	/* set segment pointers (assuming triangular facets in 3D) */  		
	void SetEdgePointers(const double* coords, int edge, const double** ppA,
		const double** ppB, const double** ppC) const;
	
	/* initialize new front node */
	FrontNodeT* NewFrontNode(const FrontSegmentT& seg1,
		const FrontSegmentT& seg2) const;

	/* resent front node */
	void ResetFrontNode(const FrontSegmentT& seg1,
		const FrontSegmentT& seg2, FrontNodeT& node) const;

	/* check front spacing - front segments shouldn't be
	 * more than fda*fda_s (sampling distance) apart */
	void CheckFrontSpacing(void);
	
private:

	/* minimum acceptable kink angle */
	double fMinAngle;

	/* crack front line segments */
	AutoArrayT<FrontSegmentT*> fFrontLines;
};

#endif /* _FRONT_3D_T_H_ */
