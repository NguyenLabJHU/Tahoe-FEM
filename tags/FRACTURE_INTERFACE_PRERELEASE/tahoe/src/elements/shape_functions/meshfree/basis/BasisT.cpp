/* $Id: BasisT.cpp,v 1.2 2002-07-02 19:57:05 cjkimme Exp $ */
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
	fArray2DGroup1(0, 0)
{
	fArray2DGroup1.Register(fP);
	
	for (int i = 0; i < fDP.Length(); i++)
		fArray2DGroup1.Register(fDP[i]);

	for (int j = 0; j < fDDP.Length(); j++)
		fArray2DGroup1.Register(fDDP[j]);
}

/***********************************************************************
* Protected
***********************************************************************/

/* dimension work space */
void BasisT::Dimension(int num_nodes)
{
	fArray2DGroup1.Dimension(BasisDimension(), num_nodes);
}
