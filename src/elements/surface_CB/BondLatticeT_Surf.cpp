/* $Id: BondLatticeT_Surf.cpp,v 1.1 2005-06-29 22:38:16 hspark Exp $ */
/* created: paklein (01/07/1997) */
#include "BondLatticeT_Surf.h"
#include <math.h>

using namespace Tahoe;

/* constructor */
BondLatticeT_Surf::BondLatticeT_Surf(void) { }

/* destructor */
BondLatticeT_Surf::~BondLatticeT_Surf(void) { }

/* initialize bond table */
void BondLatticeT_Surf::Initialize(const dMatrixT* Q)
{
	const char caller[] = "BondLatticeT_Surf::Initialize";

	/* pure virtual */
	LoadBondTable();

	/* dimension work space */
	int nsd = fBonds.MinorDim();
	fBondDp.Dimension(nsd);
//	fLatDimMatrix.Dimension(nsd);
	fStrain.Dimension(nsd);
	fStretch.Dimension(nsd);

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
int BondLatticeT_Surf::NumberOfLatticeDim(void) const { return fNumLatticeDim; }
int BondLatticeT_Surf::NumberOfSpatialDim(void) const { return fNumSpatialDim; }
int BondLatticeT_Surf::NumberOfBonds(void) const { return fNumBonds; }
#endif

/*
* Compute deformed bond lengths from the given Green strain
* tensor
*/
void BondLatticeT_Surf::ComputeDeformedLengths(const dSymMatrixT& strain)
{
	/* spatial vs. lattice dimension translation */
	if (fStrain.Rows() == 3 && strain.Rows() == 2)
		fStrain.ExpandFrom2D(strain);
	else
		fStrain = strain;

	/* compute stretch tensor */
	fStretch.SetToScaled(2.0, fStrain);
	fStretch.PlusIdentity(1.0);

	/* loop over all bonds */
	for (int bond = 0; bond < fBonds.MajorDim(); bond++)
	{
		/* get bond vector */
		fBonds.RowAlias(bond, fBondSh);
		
		/* using symmetry in C */
		fStretch.Multx(fBondSh, fBondDp);
		
		/* deformed length */
		fDefLength[bond] = sqrt(dArrayT::Dot(fBondSh, fBondDp));
	}
}
