/* $Id: BondLatticeT.cpp,v 1.6 2004-07-15 08:26:42 paklein Exp $ */
/* created: paklein (01/07/1997) */
#include "BondLatticeT.h"
#include <math.h>

using namespace Tahoe;

/* constructor */
BondLatticeT::BondLatticeT(void) { }

/* destructor */
BondLatticeT::~BondLatticeT(void) { }

/* initialize bond table */
void BondLatticeT::Initialize(const dMatrixT* Q)
{
	const char caller[] = "BondLatticeT::Initialize";

	/* pure virtual */
	LoadBondTable();

	/* dimension work space */
	int nsd = fBonds.MinorDim();
	fBondDp.Dimension(nsd);
	fLatDimMatrix.Dimension(nsd);
	fStrain.Dimension(nsd);

	/* transformation */
	if (Q) 
	{
		/* copy */
		fQ = *Q;

		/* dimension check */
		if (fQ.Rows() != fQ.Cols()) ExceptionT::GeneralFail(caller);
		if (fQ.Rows() != nsd)
			ExceptionT::SizeMismatch(caller, "Q must be %dD not %dD", nsd, fQ.Rows());
	}

	/* transform bonds */
	if ( fQ.Rows() > 0 )
		for (int bond = 0; bond < fBonds.MajorDim(); bond++)
		{
			/* get bond vector */
			fBonds.RowAlias(bond, fBondSh);
			
			/* temp */
			fBondDp = fBondSh;
		
			/* transform */
			fQ.MultTx(fBondDp, fBondSh);		
		}
}

#if 0
/* accessors */
int BondLatticeT::NumberOfLatticeDim(void) const { return fNumLatticeDim; }
int BondLatticeT::NumberOfSpatialDim(void) const { return fNumSpatialDim; }
int BondLatticeT::NumberOfBonds(void) const { return fNumBonds; }
#endif

/*
* Compute deformed bond lengths from the given Green strain
* tensor
*/
void BondLatticeT::ComputeDeformedLengths(const dSymMatrixT& strain)
{
	/* spatial vs. lattice dimension translation */
	if (fStrain.Rows() == 3 && strain.Rows() == 2)
		fStrain.ExpandFrom2D(strain);
	else
		fStrain = strain;
	
	/* compute stretch tensor */
	dMatrixT& stretch = fLatDimMatrix;
	fStrain.ToMatrix(stretch);
	stretch *= 2.0;
	stretch.PlusIdentity(1.0);

	/* loop over all bonds */
	for (int bond = 0; bond < fBonds.MajorDim(); bond++)
	{
		/* get bond vector */
		fBonds.RowAlias(bond, fBondSh);
		
		/* using symmetry in C */
		stretch.MultTx(fBondSh, fBondDp);
		
		/* deformed length */
		fDefLength[bond] = sqrt(dArrayT::Dot(fBondSh, fBondDp));
	}
}
