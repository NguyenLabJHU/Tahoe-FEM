/* $Id: FrontT.h,v 1.4 2002-09-10 13:42:44 paklein Exp $ */
/* created: paklein (02/11/2000)                                          */
/* base class for crack surface objects                                   */

#ifndef _FRONT_T_H_
#define _FRONT_T_H_

/* direct members */
#include "AutoArrayT.h"
#include "nVariArray2DT.h"
#include "dArray2DT.h"

namespace Tahoe {

/* forward declarations */
class iArrayT;
class FrontNodeT;

class FrontT
{
public:

	/* constructor */
	FrontT(int nsd, int num_facet_nodes, double cone, double da, double da_s,
		int num_pts);

	/* destructor */
	virtual ~FrontT(void);
	
	/* initialize the front with the given geometry */
	virtual void Initialize(const dArray2DT& facet_coords, const iArrayT& fr_facets,
		const iArrayT& fr_edges) = 0;
	
	/* front nodes data */
	const AutoArrayT<FrontNodeT*>& FrontNodes(void) const;
		
	/* extend front at the specified points along the given direction - returns new
	 * cutting facets */
	virtual const dArray2DT& NewFacets(const ArrayT<int>& extend_pts,
		const ArrayT<int>& extend_dir) = 0;

	/* write front nodes data to output */
	virtual void Write(ostream& out) const;
	
protected:

	/* sampling parameters */
	double fcone;
	double fda;
	double fda_s; // fraction of fda at which sampling occurs
	int    fnum_pts;

	/* crack front nodes */
	AutoArrayT<FrontNodeT*> fFrontNodes;

	/* return value */
	dArray2DT fNewFacets;
	nVariArray2DT<double> fNewFacetMan;	
};

/* inlines */
inline const AutoArrayT<FrontNodeT*>& FrontT::FrontNodes(void) const
{
	return fFrontNodes;
}

} // namespace Tahoe 
#endif /* _FRONT_T_H_ */
