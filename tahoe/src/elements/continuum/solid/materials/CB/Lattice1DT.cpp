/* $Id: Lattice1DT.cpp,v 1.1.2.1 2004-04-15 21:10:19 paklein Exp $ */
#include "Lattice1DT.h"

using namespace Tahoe;

/* constructor */
Lattice1DT::Lattice1DT(int nshells):
	CBLatticeT(1, 1, nshells),
	fNumShells(nshells)
{

}

/* initialize bond table values */
void Lattice1DT::LoadBondTable(void)
{
  	fBondCounts = 1;
  	fDefLength = 0.0; 
  
  	if (fBonds.MajorDim() != fNumBonds ||
     	fBonds.MinorDim() != 1) ExceptionT::GeneralFail("Lattice1DT::LoadBondTable");

	for (int i = 0; i < fNumShells; i++)
		fBonds(i,0) = double(i+1);
}
