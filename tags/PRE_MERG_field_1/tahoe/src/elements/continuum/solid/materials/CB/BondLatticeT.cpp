/* $Id: BondLatticeT.cpp,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (01/07/1997)                                          */
/* BondLatticeT.cpp                                                       */

#include "BondLatticeT.h"
#include <math.h>

/* Constructor */
BondLatticeT::BondLatticeT(int numlatticedim, int numspatialdim,
	int numbonds): fIsInitialized(0), fNumLatticeDim(numlatticedim),
	fNumSpatialDim(numspatialdim), fNumBonds(numbonds),
	fBondCounts(fNumBonds), fBonds(fNumBonds, fNumLatticeDim),
	fDefLength(fNumBonds), fQ(0), fBondDp(fNumLatticeDim),
	fLatDimMatrix(fNumLatticeDim), fStrain(fNumLatticeDim)
{
	/* dimension checks */
	if (fNumLatticeDim == 3)
		if (fNumSpatialDim != 2 && fNumSpatialDim != 3) throw eGeneralFail;
	else if (fNumLatticeDim == 2)
		if (fNumSpatialDim != 2) throw eGeneralFail;
	else
		throw eGeneralFail;

	if (fNumBonds < 1) throw eGeneralFail;
}

/*
* The Q matrix passed into this constructor is used to rotate the
* bond vectors into the orientation prescribed by Q.  No check is
* performed on the orthogonality of Q, only its dimensions.  Q is
* deep copied.  Q is defined as:
*
*			Q = d x_global / d x_natural
*
* So that the vectors are transformed by:
*
*			r_new = Q.r_natural
*
*/
BondLatticeT::BondLatticeT(const dMatrixT& Q, int numspatialdim,
	int numbonds): fIsInitialized(0), fNumLatticeDim(Q.Rows()),
	fNumSpatialDim(numspatialdim), fNumBonds(numbonds),
	fBondCounts(fNumBonds), fBonds(fNumBonds, fNumLatticeDim),
	fDefLength(fNumBonds), fQ(fNumLatticeDim), fBondDp(fNumLatticeDim),
	fLatDimMatrix(fNumLatticeDim), fStrain(fNumLatticeDim)
{
	/* check transformation matrix dimensions */
	if (Q.Rows() != Q.Cols()) throw eGeneralFail;
	
	/* deep copy */
	fQ = Q;

	/* dimension checks */
	if (fNumLatticeDim == 3)
		if (fNumSpatialDim != 2 && fNumSpatialDim != 3) throw eGeneralFail;
	else if (fNumLatticeDim == 2)
		if (fNumSpatialDim != 2) throw eGeneralFail;
	else
		throw eGeneralFail;

	if (fNumBonds < 1) throw eGeneralFail;
}
	
/* Destructor */
BondLatticeT::~BondLatticeT(void)
{

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
* References to the lattice lists
*/
const iArrayT& BondLatticeT::BondCounts(void) const { return fBondCounts; }
const dArrayT& BondLatticeT::DeformedLengths(void) const { return fDefLength; }

/*
* Accessors
*/
int BondLatticeT::NumberOfLatticeDim(void) const { return fNumLatticeDim; }
int BondLatticeT::NumberOfSpatialDim(void) const { return fNumSpatialDim; }
int BondLatticeT::NumberOfBonds(void) const { return fNumBonds; }

/**********************************************************************
* Protected
**********************************************************************/

/*
* Compute deformed bond lengths from the given Green strain
* tensor
*/
void BondLatticeT::ComputeDeformedLengths(const dSymMatrixT& strain)
{
	/* check fBonds has been set */
	if (!fIsInitialized) throw eGeneralFail;
	
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
