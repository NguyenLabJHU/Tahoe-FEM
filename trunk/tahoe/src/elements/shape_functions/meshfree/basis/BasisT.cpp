/* $Id: BasisT.cpp,v 1.5 2004-11-03 16:09:48 raregue Exp $ */
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
	/*  DDDP is a [nsd]x[nsd]x[nsd] or [nsd]x[nsd*nsd] matrix. 
	  	using symmetry it reduces to [nsd]x[nstr]
	  	only the first three (3D) or two (2D) columns (contribution from 
	  	diagonal terms) of the [nsd]x[nstr] matrix are needed for calculation of
	  	the Laplacian of the strain tensor
	  	DDDw, thus, becomes a [nsd]x[nsd] unsymmetric matrix 
	*/
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
