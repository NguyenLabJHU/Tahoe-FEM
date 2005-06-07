/* $Id: Lattice1DT.cpp,v 1.3 2004-07-15 08:26:42 paklein Exp $ */
#include "Lattice1DT.h"

using namespace Tahoe;

/* constructor */
Lattice1DT::Lattice1DT(int nshells):
	fNumShells(nshells)
{

}

/* initialize bond table values */
void Lattice1DT::LoadBondTable(void)
{
	/* dimension work space */
	int num_bonds = fNumShells;
	fBondCounts.Dimension(num_bonds);
	fDefLength.Dimension(num_bonds);
	fBonds.Dimension(num_bonds, 1);

	/* initialize */
  	fBondCounts = 1;
  	fDefLength = 0.0; 
  
	for (int i = 0; i < fNumShells; i++)
		fBonds(i,0) = double(i+1);
}
