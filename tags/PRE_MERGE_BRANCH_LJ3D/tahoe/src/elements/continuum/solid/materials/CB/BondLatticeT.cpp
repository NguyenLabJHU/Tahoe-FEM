/* $Id: BondLatticeT.cpp,v 1.5 2004-06-26 05:58:47 paklein Exp $ */
/* created: paklein (01/07/1997) */
#include "BondLatticeT.h"
#include <math.h>

using namespace Tahoe;

/* constructor */
BondLatticeT::BondLatticeT(int numlatticedim, int numspatialdim,
	int numbonds): fIsInitialized(0), fNumLatticeDim(numlatticedim),
	fNumSpatialDim(numspatialdim), fNumBonds(numbonds),
	fBondCounts(fNumBonds), fBonds(fNumBonds, fNumLatticeDim),
	fDefLength(fNumBonds), fQ(0), fBondDp(fNumLatticeDim),
	fLatDimMatrix(fNumLatticeDim), fStrain(fNumLatticeDim)
{
	const char caller[] = "BondLatticeT::BondLatticeT";

	/* dimension checks */
	if (fNumLatticeDim == 3)
		if (fNumSpatialDim != 2 && fNumSpatialDim != 3) ExceptionT::GeneralFail(caller);
	else if (fNumLatticeDim == 2)
		if (fNumSpatialDim != 2) ExceptionT::GeneralFail(caller);
	else
		ExceptionT::GeneralFail(caller);

	if (fNumBonds < 1) ExceptionT::GeneralFail(caller);
	
	/* initialize values */
	fBondCounts = 0;
	fDefLength = 1.0;
}

BondLatticeT::BondLatticeT(const dMatrixT& Q, int numspatialdim,
	int numbonds): fIsInitialized(0), fNumLatticeDim(Q.Rows()),
	fNumSpatialDim(numspatialdim), fNumBonds(numbonds),
	fBondCounts(fNumBonds), fBonds(fNumBonds, fNumLatticeDim),
	fDefLength(fNumBonds), fQ(fNumLatticeDim), fBondDp(fNumLatticeDim),
	fLatDimMatrix(fNumLatticeDim), fStrain(fNumLatticeDim)
{
	const char caller[] = "BondLatticeT::BondLatticeT";

	/* check transformation matrix dimensions */
	if (Q.Rows() != Q.Cols()) ExceptionT::GeneralFail(caller);
	
	/* deep copy */
	fQ = Q;

	/* dimension checks */
	if (fNumLatticeDim == 3)
		if (fNumSpatialDim != 2 && fNumSpatialDim != 3) ExceptionT::GeneralFail(caller);
	else if (fNumLatticeDim == 2)
		if (fNumSpatialDim != 2) ExceptionT::GeneralFail(caller);
	else
		ExceptionT::GeneralFail(caller);

	if (fNumBonds < 1) ExceptionT::GeneralFail(caller);
}

/* initialize bond table */
void BondLatticeT::Initialize(void)
{
	/* pure virtual */
	LoadBondTable();

	/* transform bonds */
	if ( fQ.Rows() > 0 )
		for (int bond = 0; bond < fNumBonds; bond++)
		{
			/* get bond vector */
			fBonds.RowAlias(bond, fBondSh);
			
			/* temp */
			fBondDp = fBondSh;
		
			/* transform */
			fQ.MultTx(fBondDp, fBondSh);		
		}

	/* set flag */
	fIsInitialized = 1;
}

/*
* Compute deformed bond lengths from the given Green strain
* tensor
*/
void BondLatticeT::ComputeDeformedLengths(const dSymMatrixT& strain)
{
	/* check fBonds has been set */
	if (!fIsInitialized) throw ExceptionT::kGeneralFail;
	
	/* spatial vs. lattice dimension translation */
	if (fNumLatticeDim == 3 && fNumSpatialDim == 2)
		fStrain.ExpandFrom2D(strain);
	else
		fStrain = strain;
	
	/* compute stretch tensor */
	dMatrixT& fStretch = fLatDimMatrix;
	fStrain.ToMatrix(fStretch);
	fStretch *= 2.0;
	fStretch.PlusIdentity(1.0);

	/* loop over all bonds */
	for (int bond = 0; bond < fNumBonds; bond++)
	{
		/* get bond vector */
		fBonds.RowAlias(bond, fBondSh);
		
		/* using symmetry in C */
		fStretch.MultTx(fBondSh, fBondDp);
		
		/* deformed length */
		fDefLength[bond] = sqrt(dArrayT::Dot(fBondSh, fBondDp));
	}
}
