/* $Id: FrontT.cpp,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: paklein (02/11/2000)                                          */

#include "FrontT.h"

#include "FrontNodeT.h"

/* constants */
const double     Pi = acos(-1.0);
const int kHeadRoom = 20;

/* constructor */
FrontT::FrontT(int nsd, int num_facet_nodes, double cone, double da,
	double da_s, int num_pts):
	fcone(cone),
	fda(da),
	fda_s(da_s),
	fnum_pts(num_pts),
	fFrontNodes(kHeadRoom, false),
	fNewFacetMan(20, fNewFacets, nsd*num_facet_nodes)
{
	/* checks */
	if (fcone < 0.0 || fcone > 180.0) throw eBadInputValue;
	if (fda < 0.0) throw eBadInputValue;
	if (fnum_pts < 0) throw eBadInputValue;
}

/* destructor */
FrontT::~FrontT(void)
{
	/* free front nodes data */
	for (int j = 0; j < fFrontNodes.Length(); j++)
		delete fFrontNodes[j];

	/* empty lists */
	fFrontNodes.Allocate(0);
}


/* write front nodes data to output */
void FrontT::Write(ostream& out) const
{
	for (int i = 0; i < fFrontNodes.Length(); i++)
		fFrontNodes[i]->Write(out);
}
