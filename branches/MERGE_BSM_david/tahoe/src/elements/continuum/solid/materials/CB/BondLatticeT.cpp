/* $Id: BondLatticeT.cpp,v 1.8 2006-07-03 20:19:32 hspark Exp $ */
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
	{
		for (int bond = 0; bond < fBonds.MajorDim(); bond++)
		{
			/* get bond vector */
			fBonds.RowAlias(bond, fBondSh);
			
			/* temp */
			fBondDp = fBondSh;
		
			/* transform */
			fQ.MultTx(fBondDp, fBondSh);		
		}
		/* transform bulk bonds for surface CB stuff */
		for (int bond = 0; bond < fBulkBonds.MajorDim(); bond++)
		{
			/* get bond vector */
			fBulkBonds.RowAlias(bond, fBondShB);
			
			/* temp */
			fBondDpB = fBondShB;
		
			/* transform */
			fQ.MultTx(fBondDpB, fBondShB);		
		}
		/* transform surface 1 bonds for surface CB stuff */
		for (int bond = 0; bond < fSurf1Bonds.MajorDim(); bond++)
		{
			/* get bond vector */
			fSurf1Bonds.RowAlias(bond, fBondShS1);
			
			/* temp */
			fBondDpS1 = fBondShS1;
		
			/* transform */
			fQ.MultTx(fBondDpS1, fBondShS1);		
		}
		/* transform surface 2 bonds for surface CB stuff */
		for (int bond = 0; bond < fSurf2Bonds.MajorDim(); bond++)
		{
			/* get bond vector */
			fSurf2Bonds.RowAlias(bond, fBondShS2);
			
			/* temp */
			fBondDpS2 = fBondShS2;
		
			/* transform */
			fQ.MultTx(fBondDpS2, fBondShS2);		
		}
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

/* BELOW ARE NEW SURFACE CB SPECIFIC FUNCTIONS */
/* Compute deformed lengths for a representative bulk atom */
void BondLatticeT::ComputeDeformedBulkBonds(const dSymMatrixT& strain)
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
	for (int bond = 0; bond < fBulkBonds.MajorDim(); bond++)
	{
		/* get bond vector */
		fBulkBonds.RowAlias(bond, fBondShB);
		
		/* using symmetry in C */
		fStretch.Multx(fBondShB, fBondDpB);
		
		/* deformed length */
		fDefBulk[bond] = sqrt(dArrayT::Dot(fBondShB, fBondDpB));
	}
}

/* Compute deformed lengths for a representative surface atom */
void BondLatticeT::ComputeDeformedSurf1Bonds(const dSymMatrixT& strain)
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
	for (int bond = 0; bond < fSurf1Bonds.MajorDim(); bond++)
	{
		/* get bond vector */
		fSurf1Bonds.RowAlias(bond, fBondShS1);
		
		/* using symmetry in C */
		fStretch.Multx(fBondShS1, fBondDpS1);
		
		/* deformed length */
		fDefSurf1[bond] = sqrt(dArrayT::Dot(fBondShS1, fBondDpS1));
	}
}

/* Compute deformed lengths for a representative atom 1 layer into the bulk */
void BondLatticeT::ComputeDeformedSurf2Bonds(const dSymMatrixT& strain)
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
	for (int bond = 0; bond < fSurf2Bonds.MajorDim(); bond++)
	{
		/* get bond vector */
		fSurf2Bonds.RowAlias(bond, fBondShS2);
		
		/* using symmetry in C */
		fStretch.Multx(fBondShS2, fBondDpS2);
		
		/* deformed length */
		fDefSurf2[bond] = sqrt(dArrayT::Dot(fBondShS2, fBondDpS2));
	}
}
