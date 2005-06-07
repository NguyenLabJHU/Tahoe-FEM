/* $Id: HexLattice2DT.cpp,v 1.3 2004-07-15 08:26:42 paklein Exp $ */
#include "HexLattice2DT.h"

using namespace Tahoe;

/* number of atoms per shell */
const int atoms_per_shell[] = {3, 3, 3, 6, 3};
const int atoms_in_shells[] = {3, 6, 9, 15, 18};
static int AtomsInShells(int nshells) {
	if (nshells < 0 || nshells > 5) ExceptionT::OutOfRange();
	return atoms_in_shells[nshells-1];
};
const double sqrt3 = sqrt(3.0);

/* constructor */
HexLattice2DT::HexLattice2DT(int nshells):
	fNumShells(nshells)
{

}

/* initialize bond table values */
void HexLattice2DT::LoadBondTable(void)
{
	/* dimension work space */
	int num_bonds = AtomsInShells(fNumShells);
	fBondCounts.Dimension(num_bonds);
	fDefLength.Dimension(num_bonds);
	fBonds.Dimension(num_bonds, 2);

	/* initialize */
  	fBondCounts = 1;
  	fDefLength = 0.0; 
  
  	double bonddata1[3*2] =
  		{ 1.0, 0.0,
  		  0.5, sqrt3/2.0,
  		 -0.5, sqrt3/2.0};

  	double bonddata2[3*2] =
  		{ 0.0, sqrt3,
  		  1.5, sqrt3/2.0,
  		 -1.5, sqrt3/2.0};

  	double bonddata3[3*2] =
  		{ 2.0, 0.0,
  		  1.0, sqrt3,
  		 -1.0, sqrt3};

  	double bonddata4[6*2] =
  		{ 2.5, sqrt3/2.0,
  		  2.0, sqrt3,
  		  0.5, 3.0*sqrt3/2.0,
  		 -0.5, 3.0*sqrt3/2.0,
  		 -2.0, sqrt3,
  		 -2.5, sqrt3/2.0};

  	double bonddata5[3*2] =
  		{ 3.0, 0.0,
  		  1.5, 3.0*sqrt3/2.0,
  		 -1.5, 3.0*sqrt3/2.0};

	double* shells[5];
	shells[0] = bonddata1;
	shells[1] = bonddata2;
	shells[2] = bonddata3;
	shells[3] = bonddata4;
	shells[4] = bonddata5;

	int bond = 0;
	for (int i = 0; i < fNumShells; i++)
	{
		dArray2DT bonds(atoms_per_shell[i], 2, shells[i]);
		for (int j = 0; j < bonds.MajorDim(); j++)
		{
			fBonds(bond,0) = bonds(j,0);
			fBonds(bond,1) = bonds(j,1);
			bond++;
		}
	}
}
