/* $Id: BasisT.cpp,v 1.1 2004-08-14 00:03:35 raregue Exp $ */
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
	//fDDDP(fNumSD*dSymMatrixT::NumValues(fNumSD)), // kyonten
	fDDDP(fNumSD*fNumSD), // kyonten
	// DDDp(nsd*nsd) is a special case. In 3D case out of
	// 27 components, only 9 are needed
	fArray2DGroup1(0, 0)
{
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
