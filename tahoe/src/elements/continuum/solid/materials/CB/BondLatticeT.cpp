/* $Id: BondLatticeT.cpp,v 1.9 2009-06-06 17:28:49 hspark Exp $ */
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
	fBondDpS.Dimension(nsd);
	fBondDpE.Dimension(nsd);
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
		/* TRANSFORM BONDS FOR EDGE CB STUFF */
		for (int bond = 0; bond < fAllBulkBonds.MajorDim(); bond++)
		{
			/* get bond vector */
			fAllBulkBonds.RowAlias(bond, fBondShBE);
			
			/* temp */
			fBondDpBE = fBondShBE;
		
			/* transform */
			fQ.MultTx(fBondDpBE, fBondShBE);		
		}
		/* Edge 1 atoms */
		for (int bond = 0; bond < fEdgeBonds1.MajorDim(); bond++)
		{
			/* get bond vector */
			fEdgeBonds1.RowAlias(bond, fBondShE1);
			
			/* temp */
			fBondDpE1 = fBondShE1;
		
			/* transform */
			fQ.MultTx(fBondDpE1, fBondShE1);		
		}
		/* Edge 2 atoms */
		for (int bond = 0; bond < fEdgeBonds2.MajorDim(); bond++)
		{
			/* get bond vector */
			fEdgeBonds2.RowAlias(bond, fBondShE2);
			
			/* temp */
			fBondDpE2 = fBondShE2;
		
			/* transform */
			fQ.MultTx(fBondDpE2, fBondShE2);		
		}
		/* Edge 3 atoms */
		for (int bond = 0; bond < fEdgeBonds3.MajorDim(); bond++)
		{
			/* get bond vector */
			fEdgeBonds3.RowAlias(bond, fBondShE3);
			
			/* temp */
			fBondDpE3 = fBondShE3;
		
			/* transform */
			fQ.MultTx(fBondDpE3, fBondShE3);		
		}
		/* Edge 4 atoms */
		for (int bond = 0; bond < fEdgeBonds4.MajorDim(); bond++)
		{
			/* get bond vector */
			fEdgeBonds4.RowAlias(bond, fBondShE4);
			
			/* temp */
			fBondDpE4 = fBondShE4;
		
			/* transform */
			fQ.MultTx(fBondDpE4, fBondShE4);		
		}
		/* Edge 5 atoms */
		for (int bond = 0; bond < fEdgeBonds5.MajorDim(); bond++)
		{
			/* get bond vector */
			fEdgeBonds5.RowAlias(bond, fBondShE5);
			
			/* temp */
			fBondDpE5 = fBondShE5;
		
			/* transform */
			fQ.MultTx(fBondDpE5, fBondShE5);		
		}
		/* Edge 6 atoms */
		for (int bond = 0; bond < fEdgeBonds6.MajorDim(); bond++)
		{
			/* get bond vector */
			fEdgeBonds6.RowAlias(bond, fBondShE6);
			
			/* temp */
			fBondDpE6 = fBondShE6;
		
			/* transform */
			fQ.MultTx(fBondDpE6, fBondShE6);		
		}
		/* Edge 7 atoms */
		for (int bond = 0; bond < fEdgeBonds7.MajorDim(); bond++)
		{
			/* get bond vector */
			fEdgeBonds7.RowAlias(bond, fBondShE7);
			
			/* temp */
			fBondDpE7 = fBondShE7;
		
			/* transform */
			fQ.MultTx(fBondDpE7, fBondShE7);		
		}
		/* Edge 8 atoms */
		for (int bond = 0; bond < fEdgeBonds8.MajorDim(); bond++)
		{
			/* get bond vector */
			fEdgeBonds8.RowAlias(bond, fBondShE8);
			
			/* temp */
			fBondDpE8 = fBondShE8;
		
			/* transform */
			fQ.MultTx(fBondDpE8, fBondShE8);		
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

/* Compute deformed surface bond lengths */
void BondLatticeT::ComputeDeformedSurfaceLengths(const dSymMatrixT& strain)
{
	/* spatial vs. lattice dimension translation */
	if (fStrain.Rows() == 3 && strain.Rows() == 2)
		fStrain.ExpandFrom2D(strain);
	else
		fStrain = strain;

	/* compute stretch tensor */
	fStretch.SetToScaled(2.0, fStrain);
	fStretch.PlusIdentity(1.0);
//	fStretch.Identity();

	/* loop over all bonds */
	for (int bond = 0; bond < fBonds.MajorDim(); bond++)
	{
		/* get bond vector */
		fBonds.RowAlias(bond, fBondShS);

		/* using symmetry in C */
		fStretch.Multx(fBondShS, fBondDpS);

		/* deformed length */
		fDefLength[bond] = sqrt(dArrayT::Dot(fBondShS, fBondDpS));
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
//	fStretch.Identity();

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
//	fStretch.Identity();

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
//	fStretch.Identity();

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

/* Compute deformed EDGE bond lengths */
void BondLatticeT::ComputeDeformedEdgeLengths(const dSymMatrixT& strain)
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
		fBonds.RowAlias(bond, fBondShE);

		/* using symmetry in C */
		fStretch.Multx(fBondShE, fBondDpE);

		/* deformed length */
		fDefLength[bond] = sqrt(dArrayT::Dot(fBondShE, fBondDpE));
	}
}

void BondLatticeT::ComputeDeformedEdgeBulkBonds(const dSymMatrixT& strain)
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
	for (int bond = 0; bond < fAllBulkBonds.MajorDim(); bond++)
	{
		/* get bond vector */
		fAllBulkBonds.RowAlias(bond, fBondShBE);
		
		/* using symmetry in C */
		fStretch.Multx(fBondShBE, fBondDpBE);
		
		/* deformed length */
		fDefBulkLength[bond] = sqrt(dArrayT::Dot(fBondShBE, fBondDpBE));
	}
}

/* Type-1 Edge Atoms */
void BondLatticeT::ComputeDeformedEdge1Bonds(const dSymMatrixT& strain)
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
	for (int bond = 0; bond < fEdgeBonds1.MajorDim(); bond++)
	{
		/* get bond vector */
		fEdgeBonds1.RowAlias(bond, fBondShE1);
		
		/* using symmetry in C */
		fStretch.Multx(fBondShE1, fBondDpE1);
		
		/* deformed length */
		fDefEdge1[bond] = sqrt(dArrayT::Dot(fBondShE1, fBondDpE1));
	}
}

/* Type-2 Edge Atoms */
void BondLatticeT::ComputeDeformedEdge2Bonds(const dSymMatrixT& strain)
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
	for (int bond = 0; bond < fEdgeBonds2.MajorDim(); bond++)
	{
		/* get bond vector */
		fEdgeBonds2.RowAlias(bond, fBondShE2);
		
		/* using symmetry in C */
		fStretch.Multx(fBondShE2, fBondDpE2);
		
		/* deformed length */
		fDefEdge2[bond] = sqrt(dArrayT::Dot(fBondShE2, fBondDpE2));
	}
}

/* Type-3 Edge Atoms */
void BondLatticeT::ComputeDeformedEdge3Bonds(const dSymMatrixT& strain)
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
	for (int bond = 0; bond < fEdgeBonds3.MajorDim(); bond++)
	{
		/* get bond vector */
		fEdgeBonds3.RowAlias(bond, fBondShE3);
		
		/* using symmetry in C */
		fStretch.Multx(fBondShE3, fBondDpE3);
		
		/* deformed length */
		fDefEdge3[bond] = sqrt(dArrayT::Dot(fBondShE3, fBondDpE3));
	}
}

/* Type-4 Edge Atoms */
void BondLatticeT::ComputeDeformedEdge4Bonds(const dSymMatrixT& strain)
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
	for (int bond = 0; bond < fEdgeBonds4.MajorDim(); bond++)
	{
		/* get bond vector */
		fEdgeBonds4.RowAlias(bond, fBondShE4);
		
		/* using symmetry in C */
		fStretch.Multx(fBondShE4, fBondDpE4);
		
		/* deformed length */
		fDefEdge4[bond] = sqrt(dArrayT::Dot(fBondShE4, fBondDpE4));
	}
}

/* Type-5 Edge Atoms */
void BondLatticeT::ComputeDeformedEdge5Bonds(const dSymMatrixT& strain)
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
	for (int bond = 0; bond < fEdgeBonds5.MajorDim(); bond++)
	{
		/* get bond vector */
		fEdgeBonds5.RowAlias(bond, fBondShE5);
		
		/* using symmetry in C */
		fStretch.Multx(fBondShE5, fBondDpE5);
		
		/* deformed length */
		fDefEdge5[bond] = sqrt(dArrayT::Dot(fBondShE5, fBondDpE5));
	}
}

/* Type-6 Edge Atoms */
void BondLatticeT::ComputeDeformedEdge6Bonds(const dSymMatrixT& strain)
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
	for (int bond = 0; bond < fEdgeBonds6.MajorDim(); bond++)
	{
		/* get bond vector */
		fEdgeBonds6.RowAlias(bond, fBondShE6);
		
		/* using symmetry in C */
		fStretch.Multx(fBondShE6, fBondDpE6);
		
		/* deformed length */
		fDefEdge6[bond] = sqrt(dArrayT::Dot(fBondShE6, fBondDpE6));
	}
}

/* Type-7 Edge Atoms */
void BondLatticeT::ComputeDeformedEdge7Bonds(const dSymMatrixT& strain)
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
	for (int bond = 0; bond < fEdgeBonds7.MajorDim(); bond++)
	{
		/* get bond vector */
		fEdgeBonds7.RowAlias(bond, fBondShE7);
		
		/* using symmetry in C */
		fStretch.Multx(fBondShE7, fBondDpE7);
		
		/* deformed length */
		fDefEdge7[bond] = sqrt(dArrayT::Dot(fBondShE7, fBondDpE7));
	}
}

/* Type-8 Edge Atoms */
void BondLatticeT::ComputeDeformedEdge8Bonds(const dSymMatrixT& strain)
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
	for (int bond = 0; bond < fEdgeBonds8.MajorDim(); bond++)
	{
		/* get bond vector */
		fEdgeBonds8.RowAlias(bond, fBondShE8);
		
		/* using symmetry in C */
		fStretch.Multx(fBondShE8, fBondDpE8);
		
		/* deformed length */
		fDefEdge8[bond] = sqrt(dArrayT::Dot(fBondShE8, fBondDpE8));
	}
}