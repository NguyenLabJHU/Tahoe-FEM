/* $Id: FCCLatticeT.cpp,v 1.2 2003-03-31 23:14:38 paklein Exp $ */
#include "FCCLatticeT.h"

using namespace Tahoe;

/* number of atoms per shell */
const int atoms_per_shell[] = {6, 3, 12, 6, 12};
const int atoms_in_shells[] = {6, 9, 21, 27, 39};
static int AtomsInShells(int nshells) {
	if (nshells < 0 || nshells > 5) ExceptionT::OutOfRange();
	return atoms_in_shells[nshells-1];
};
const double sqrt2 = sqrt(2.0);

/* constructor */
FCCLatticeT::FCCLatticeT(const dMatrixT& Q, int nshells):
	CBLatticeT(Q, 3, AtomsInShells(nshells)),
	fNumShells(nshells)
{

}

/* initialize bond table values */
void FCCLatticeT::LoadBondTable(void)
{
  	fBondCounts = 1;
  	fDefLength = 0.0; 

  	double bonddata1[6*3] =
	{ 1.0/sqrt2, 1.0/sqrt2,       0.0,
	 -1.0/sqrt2, 1.0/sqrt2,       0.0,
	  1.0/sqrt2,       0.0, 1.0/sqrt2,
	 -1.0/sqrt2,       0.0, 1.0/sqrt2,
	        0.0, 1.0/sqrt2, 1.0/sqrt2,
	        0.0,-1.0/sqrt2, 1.0/sqrt2};

  	double bonddata2[3*3] =
	{sqrt2,   0.0,   0.0,
	   0.0, sqrt2,   0.0,
	   0.0,   0.0, sqrt2};

  	double bonddata3[12*3] =
	{     sqrt2, 1.0/sqrt2, 1.0/sqrt2,
      1.0/sqrt2,     sqrt2, 1.0/sqrt2,
      1.0/sqrt2, 1.0/sqrt2,     sqrt2,
         -sqrt2, 1.0/sqrt2, 1.0/sqrt2,
     -1.0/sqrt2,     sqrt2, 1.0/sqrt2,
     -1.0/sqrt2, 1.0/sqrt2,     sqrt2,
          sqrt2,-1.0/sqrt2, 1.0/sqrt2,
      1.0/sqrt2,    -sqrt2, 1.0/sqrt2,
      1.0/sqrt2,-1.0/sqrt2,     sqrt2,
         -sqrt2,-1.0/sqrt2, 1.0/sqrt2,
     -1.0/sqrt2,    -sqrt2, 1.0/sqrt2,
     -1.0/sqrt2,-1.0/sqrt2,     sqrt2};

  	double bonddata4[6*3] =
	{ sqrt2, sqrt2,  0.0,
     -sqrt2, sqrt2,  0.0,
      sqrt2,   0.0, sqrt2,
     -sqrt2,   0.0, sqrt2,
        0.0, sqrt2, sqrt2,
        0.0,-sqrt2, sqrt2};

  	double bonddata5[12*3] =
	{1.0/sqrt2, 3.0/sqrt2,      0.0,
    -1.0/sqrt2, 3.0/sqrt2,      0.0,
     3.0/sqrt2, 1.0/sqrt2,      0.0,
     3.0/sqrt2,-1.0/sqrt2,      0.0,
     1.0/sqrt2,       0.0, 3.0/sqrt2,
    -1.0/sqrt2,       0.0, 3.0/sqrt2,
     3.0/sqrt2,       0.0, 1.0/sqrt2,
    -3.0/sqrt2,       0.0, 1.0/sqrt2,
           0.0, 3.0/sqrt2, 1.0/sqrt2,
           0.0,-3.0/sqrt2, 1.0/sqrt2,
           0.0, 1.0/sqrt2, 3.0/sqrt2,
           0.0,-1.0/sqrt2, 3.0/sqrt2};

	double* shells[5];

	shells[0] = bonddata1;
	shells[1] = bonddata2;
	shells[2] = bonddata3;
	shells[3] = bonddata4;
	shells[4] = bonddata5;

  	if (fBonds.MajorDim() != fNumBonds ||
     	fBonds.MinorDim() != 3) ExceptionT::GeneralFail();

	int bond = 0;
	for (int i = 0; i < fNumShells; i++)
	{
		dArray2DT bonds(atoms_per_shell[i], 3, shells[i]);
		for (int j = 0; j < bonds.MajorDim(); j++)
		{
			fBonds(bond,0) = bonds(j,0);
			fBonds(bond,1) = bonds(j,1);
			fBonds(bond,2) = bonds(j,2);
			bond++;
		}
	}
}
