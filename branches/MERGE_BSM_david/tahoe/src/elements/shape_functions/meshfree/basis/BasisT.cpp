/* $Id: BasisT.cpp,v 1.6 2005-12-23 03:24:00 kyonten Exp $ */
/* created: paklein (12/10/1999)                                          */
/* base class for basis functions                                         */

#include "BasisT.h"
#include "dSymMatrixT.h"

/* constructor */

using namespace Tahoe;

BasisT::BasisT(int complete, int nsd):
	fComplete(complete),
	fNumSD(nsd),
	fDP(fNumSD),
	fDDP(dSymMatrixT::NumValues(fNumSD)),
	fDDDP(fNumSD*fNumSD), // kyonten
	fArray2DGroup1(0, 0)
{
	/* fDDDP dimension in 3D */
	if(fNumSD == 3) fDDDP.Dimension(fNumSD*fNumSD+1);
	
	fArray2DGroup1.Register(fP);
	
	for (int i = 0; i < fDP.Length(); i++)
		fArray2DGroup1.Register(fDP[i]);

	for (int j = 0; j < fDDP.Length(); j++)
		fArray2DGroup1.Register(fDDP[j]);
		
	for (int k = 0; k < fDDDP.Length(); k++) // kyonten
		fArray2DGroup1.Register(fDDDP[k]);	
}

/***********************************************************************
* Protected
***********************************************************************/

/* dimension work space */
void BasisT::Dimension(int num_nodes)
{
	fArray2DGroup1.Dimension(BasisDimension(), num_nodes);
}
