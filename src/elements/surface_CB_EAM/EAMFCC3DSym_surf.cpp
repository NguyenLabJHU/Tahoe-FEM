/* $Id: EAMFCC3DSym_surf.cpp,v 1.4 2006-07-03 20:20:09 hspark Exp $ */
/* created: paklein (12/06/1996) */
#include "EAMFCC3DSym_surf.h"

using namespace Tahoe;

/* Bond table parameters */
const int kEAMFCC3DSurfBonds        = 78;
const int kEAMFCC3DNumBonds			= 54;
const int kEAMFCC3DSurf1Bonds       = 33;
const int kEAMFCC3DSurf2Bonds       = 45;
const int kEAMFCC3DNumLatticeDim 	=  3;
const int kEAMFCC3DNumAtomsPerCell	=  4;
const int kEAMFCC3DNumAtomsPerArea  =  2;
const double piby2 = 4.0 * atan(1.0) / 2.0;

/* constructor */
EAMFCC3DSym_surf::EAMFCC3DSym_surf(int nshells, int normal):
	EAMFCC3D_surf(nshells, normal),
	fNormalCode(normal)
{
	SetName("FCC_EAM_Cauchy-Born");
}

/**********************************************************************
 * Protected
 **********************************************************************/
	
void EAMFCC3DSym_surf::LoadBondTable(void)
{
	/* dimension work space - ARE THESE DIMENSIONS CORRECT? */
	fBondCounts.Dimension(kEAMFCC3DSurfBonds);
	fBulkCounts.Dimension(kEAMFCC3DNumBonds);
	fSurf1Counts.Dimension(kEAMFCC3DSurf1Bonds);
	fSurf2Counts.Dimension(kEAMFCC3DSurf2Bonds);
	fDefLength.Dimension(kEAMFCC3DSurfBonds);
	fDefBulk.Dimension(kEAMFCC3DNumBonds);
	fDefSurf1.Dimension(kEAMFCC3DSurf1Bonds);
	fDefSurf2.Dimension(kEAMFCC3DSurf2Bonds);
	fBonds.Dimension(kEAMFCC3DSurfBonds, kEAMFCC3DNumLatticeDim);
	fBulkBonds.Dimension(kEAMFCC3DNumBonds,3);
	fSurf1Bonds.Dimension(kEAMFCC3DSurf1Bonds,3);
	fSurf2Bonds.Dimension(kEAMFCC3DSurf2Bonds,3);
	fAtomType.Dimension(kEAMFCC3DSurfBonds);

	dArray2DT temp_bonds, temp_bonds2;
	temp_bonds.Dimension(kEAMFCC3DSurfBonds, 3);	// temporary bond table before rotation
	temp_bonds2.Dimension(kEAMFCC3DSurfBonds, 3);	// Currently have # of bonds for {100} surfaces

	/* all bonds appear once */
	fBondCounts = 1;
	fBulkCounts = 1;
	fSurf1Counts = 1;
	fSurf2Counts = 1;
	
	/* clear deformed lengths for now */
	fDefLength = 0.0;
	fDefBulk = 0.0;
	fDefSurf1 = 0.0;
	fDefSurf2 = 0.0;

	/* undeformed bond data for bulk atom with 4th neighbor interactions */
	double bulkbond[kEAMFCC3DNumBonds][kEAMFCC3DNumLatticeDim] = {
		{0, 0, -1.},
		{0, 0, 1.},
		{0, -1., 0},
		{0, 1., 0},
		{-1., 0, 0},
		{1., 0, 0},
		{-0.5, 0, -0.5},
		{-0.5, 0, 0.5},
		{0.5, 0, -0.5},
		{0.5, 0, 0.5},
		{0, -0.5, -0.5},
		{0, -0.5, 0.5},
		{0, 0.5, -0.5},
		{0, 0.5, 0.5},
		{-0.5, -0.5, 0},
		{-0.5, 0.5, 0},
		{0.5, -0.5, 0},
		{0.5, 0.5, 0},
		{-1., 0, -1.},
		{-1., 0, 1.},
		{1., 0, -1.},
		{1., 0, 1.},
		{0, -1., -1.},
		{0, -1., 1.},
		{0, 1., -1.},
		{0, 1., 1.},
		{-1., -1., 0},
		{-1., 1., 0},
		{1., -1., 0},
		{1., 1., 0},
		{-0.5, -0.5, -1.},
		{-0.5, -0.5, 1.},
		{-0.5, 0.5, -1.},
		{-0.5, 0.5, 1.},
		{0.5, -0.5, -1.},
		{0.5, -0.5, 1.},
		{0.5, 0.5, -1.},
		{0.5, 0.5, 1.},
		{-0.5, -1., -0.5},
		{-0.5, -1., 0.5},
		{-0.5, 1., -0.5},
		{-0.5, 1., 0.5},
		{0.5, -1., -0.5},
		{0.5, -1., 0.5},
		{0.5, 1., -0.5},
		{0.5, 1., 0.5},
		{-1., -0.5, -0.5},
		{-1., -0.5, 0.5},
		{-1., 0.5, -0.5},
		{-1., 0.5, 0.5},
		{1., -0.5, -0.5},
		{1., -0.5, 0.5},
		{1., 0.5, -0.5},
		{1., 0.5, 0.5}
	};

	/* Copy bond table into array */
	for (int i = 0; i < kEAMFCC3DNumBonds; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			fBulkBonds(i,j) = bulkbond[i][j];

	/* Bond table for an atom on the surface */
	double surf1bond[kEAMFCC3DSurf1Bonds][kEAMFCC3DNumLatticeDim] = {
		{0.5, 0.5, 0.0}, // Surface cluster (8 nearest neighbors)
		{0.5, -0.5, 0.0},
		{0.5, 0.0, 0.5},
		{0.5, 0.0, -0.5},
		{0.0, 0.5, 0.5},
		{0.0, -0.5, 0.5},
		{0.0, 0.5, -0.5},
		{0.0, -0.5, -0.5},
		{1.0, 0.0, 0.0}, // Surface cluster (5 2nd shell neighbors)
		{0.0, 1.0, 0.0},
		{0.0, 0.0, 1.0},
		{0.0, -1.0, 0.0},
		{0.0, 0.0, -1.0},
		{1.0, 0.5, 0.5}, // Surface cluster (12 3rd shell neighbors)
		{0.5, 1.0, 0.5},
		{0.5, 0.5, 1.0},
		{1.0, 0.5, -0.5},
		{0.5, 1.0, -0.5},
		{0.5, 0.5, -1.0},
		{1.0, -0.5, 0.5},
		{0.5, -1.0, 0.5},
		{0.5, -0.5, 1.0},
		{1.0, -0.5, -0.5},
		{0.5, -1.0, -0.5},
		{0.5, -0.5, -1.0},
		{0.0, 1.0, -1.0}, // Surface cluster (8 4th shell neighbors)
		{0.0, 1.0, 1.0},
		{1.0, 1.0, 0.0},
		{1.0, -1.0, 0.0},
		{1.0, 0.0, 1.0},
		{1.0, 0.0, -1.0},
		{0.0, -1.0, -1.0},
		{0.0, -1.0, 1.0}
	};

	/* Copy bond table into array */
	for (int i = 0; i < kEAMFCC3DSurf1Bonds; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			fSurf1Bonds(i,j) = surf1bond[i][j];

	/* Bond table for an atom 1 layer into the bulk */
	double surf2bond[kEAMFCC3DSurf2Bonds][kEAMFCC3DNumLatticeDim] = {
		{0.5, 0.5, 0.0}, // Surface cluster (8 nearest neighbors)
		{0.5, -0.5, 0.0},
		{0.5, 0.0, 0.5},
		{0.5, 0.0, -0.5},
		{0.0, 0.5, 0.5},
		{0.0, -0.5, 0.5},
		{0.0, 0.5, -0.5},
		{0.0, -0.5, -0.5},
		{-0.5, -0.5, 0.0}, // New bonds for second surface cluster begin here
		{-0.5, 0.5, 0.0},
		{-0.5, 0.0, 0.5},
		{-0.5, 0.0, -0.5},
		{1.0, 0.0, 0.0}, // Surface cluster (5 2nd shell neighbors)
		{0.0, 1.0, 0.0},
		{0.0, 0.0, 1.0},
		{0.0, -1.0, 0.0},
		{0.0, 0.0, -1.0},
		{1.0, 0.5, 0.5}, // Surface cluster (12 3rd shell neighbors)
		{0.5, 1.0, 0.5},
		{0.5, 0.5, 1.0},
		{1.0, 0.5, -0.5},
		{0.5, 1.0, -0.5},
		{0.5, 0.5, -1.0},
		{1.0, -0.5, 0.5},
		{0.5, -1.0, 0.5},
		{0.5, -0.5, 1.0},
		{1.0, -0.5, -0.5},
		{0.5, -1.0, -0.5},
		{0.5, -0.5, -1.0},
		{-0.5, 1.0, 0.5}, // New bonds for second surface cluster begin here
		{-0.5, 1.0, -0.5},
		{-0.5, 0.5, 1.0},
		{-0.5, 0.5, -1.0},
		{-0.5, -0.5, 1.0},
		{-0.5, -0.5, -1.0},
		{-0.5, -1.0, 0.5},
		{-0.5, -1.0, -0.5},
		{0.0, 1.0, -1.0}, // Surface cluster (8 4th shell neighbors)
		{0.0, 1.0, 1.0},
		{1.0, 1.0, 0.0},
		{1.0, -1.0, 0.0},
		{1.0, 0.0, 1.0},
		{1.0, 0.0, -1.0},
		{0.0, -1.0, -1.0},
		{0.0, -1.0, 1.0}
	};

	/* Copy bond table into array */
	for (int i = 0; i < kEAMFCC3DSurf2Bonds; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			fSurf2Bonds(i,j) = surf2bond[i][j];

	/* work space arrays for storing interaction types */
	iArrayT surf1nn(8);
	iArrayT surf2nn(12);
	iArrayT surf12nn(5);
	iArrayT surf22nn(5);
	iArrayT surf13nn(12);
	iArrayT surf23nn(20);
	iArrayT surf14nn(8);
	iArrayT surf24nn(8);
	int surf1n[8]={1,1,1,1,0,0,0,0};
	int surf2n[12]={5,5,5,5,4,4,4,4,3,3,3,3};
	int surf12n[5]={2,0,0,0,0};
	int surf22n[5]={5,4,4,4,4};
	int surf13n[12]={2,1,1,2,1,1,2,1,1,2,1,1};
	int surf23n[20]={5,5,5,5,5,5,5,5,5,5,5,5,3,3,3,3,3,3,3,3};
	int surf14n[8]={0,0,2,2,2,2,0,0};
	int surf24n[8]={4,4,5,5,5,5,4,4};
	surf1nn = surf1n;
	surf2nn = surf2n;
	surf12nn = surf12n;
	surf22nn = surf22n;
	surf13nn = surf13n;
	surf23nn = surf23n;
	surf14nn = surf14n;
	surf24nn = surf24n;
	fAtomType.CopyIn(0, surf1nn);
	fAtomType.CopyIn(surf1nn.Length(), surf2nn);
	fAtomType.CopyIn(surf1nn.Length()+surf2nn.Length(), surf12nn);
	fAtomType.CopyIn(surf1nn.Length()+surf2nn.Length()+surf12nn.Length(), surf22nn);
	fAtomType.CopyIn(surf1nn.Length()+surf2nn.Length()+surf12nn.Length()+surf22nn.Length(), surf13nn);
	int blah = surf1nn.Length()+surf2nn.Length()+surf12nn.Length()+surf22nn.Length()+surf13nn.Length();
	fAtomType.CopyIn(blah, surf23nn);
	fAtomType.CopyIn(blah+surf23nn.Length(), surf14nn);
	fAtomType.CopyIn(blah+surf23nn.Length()+surf14nn.Length(), surf24nn);
	
	/* New bond table for surface clusters - change dimensions! */
	/* AVOID HARD CODING NUMBER OF BONDS SPECIFIC TO {100} */
	double bonddata[kEAMFCC3DSurfBonds][kEAMFCC3DNumLatticeDim] = {
		{0.5, 0.5, 0.0}, // Surface cluster (8 nearest neighbors)
		{0.5, -0.5, 0.0},
		{0.5, 0.0, 0.5},
		{0.5, 0.0, -0.5},
		{0.0, 0.5, 0.5},
		{0.0, -0.5, 0.5},
		{0.0, 0.5, -0.5},
		{0.0, -0.5, -0.5},
		{0.5, 0.5, 0.0}, // Repeat cluster here - one atomic thickness into bulk
		{0.5, -0.5, 0.0}, // Total of 12 nearest neighbors
		{0.5, 0.0, 0.5},
		{0.5, 0.0, -0.5},
		{0.0, 0.5, 0.5},
		{0.0, -0.5, 0.5},
		{0.0, 0.5, -0.5},
		{0.0, -0.5, -0.5},
		{-0.5, -0.5, 0.0}, // New bonds for second surface cluster begin here
		{-0.5, 0.5, 0.0},
		{-0.5, 0.0, 0.5},
		{-0.5, 0.0, -0.5},
		{1.0, 0.0, 0.0}, // Surface cluster (5 2nd shell neighbors)
		{0.0, 1.0, 0.0},
		{0.0, 0.0, 1.0},
		{0.0, -1.0, 0.0},
		{0.0, 0.0, -1.0},
		{1.0, 0.0, 0.0}, // Repeat cluster here - one atomic thickness into bulk
		{0.0, 1.0, 0.0}, // Total of 5 2nd shell neighbors
		{0.0, 0.0, 1.0},
		{0.0, -1.0, 0.0},
		{0.0, 0.0, -1.0},
		{1.0, 0.5, 0.5}, // Surface cluster (12 3rd shell neighbors)
		{0.5, 1.0, 0.5},
		{0.5, 0.5, 1.0},
		{1.0, 0.5, -0.5},
		{0.5, 1.0, -0.5},
		{0.5, 0.5, -1.0},
		{1.0, -0.5, 0.5},
		{0.5, -1.0, 0.5},
		{0.5, -0.5, 1.0},
		{1.0, -0.5, -0.5},
		{0.5, -1.0, -0.5},
		{0.5, -0.5, -1.0},
		{1.0, 0.5, 0.5}, // Repeat cluster here - one atomic thickness into bulk
		{0.5, 1.0, 0.5}, // Total of 20 3rd shell neighbors
		{0.5, 0.5, 1.0},
		{1.0, 0.5, -0.5},
		{0.5, 1.0, -0.5},
		{0.5, 0.5, -1.0},
		{1.0, -0.5, 0.5},
		{0.5, -1.0, 0.5},
		{0.5, -0.5, 1.0},
		{1.0, -0.5, -0.5},
		{0.5, -1.0, -0.5},
		{0.5, -0.5, -1.0},
		{-0.5, 1.0, 0.5}, // New bonds for second surface cluster begin here
		{-0.5, 1.0, -0.5},
		{-0.5, 0.5, 1.0},
		{-0.5, 0.5, -1.0},
		{-0.5, -0.5, 1.0},
		{-0.5, -0.5, -1.0},
		{-0.5, -1.0, 0.5},
		{-0.5, -1.0, -0.5},
		{0.0, 1.0, -1.0}, // Surface cluster (8 4th shell neighbors)
		{0.0, 1.0, 1.0},
		{1.0, 1.0, 0.0},
		{1.0, -1.0, 0.0},
		{1.0, 0.0, 1.0},
		{1.0, 0.0, -1.0},
		{0.0, -1.0, -1.0},
		{0.0, -1.0, 1.0},
		{0.0, 1.0, -1.0}, // Repeat cluster here - one atomic thickness into the bulk
		{0.0, 1.0, 1.0},
		{1.0, 1.0, 0.0},
		{1.0, -1.0, 0.0},
		{1.0, 0.0, 1.0},
		{1.0, 0.0, -1.0},
		{0.0, -1.0, -1.0},
		{0.0, -1.0, 1.0}
	};

	/* Define interaction type (0-5) based on bond table */
	
	/* Rotate Bond Tables based on fNormalCode and rotation matrices */
	/* Create temporary bond table temp_bonds that combines bonddata */
	for (int i = 0; i < kEAMFCC3DSurfBonds; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			temp_bonds(i,j) = bonddata[i][j];
	
	/* Now manipulate temp_bonds */
	dMatrixT blah1(3);
	dArrayT asdf(3), prod(3);
	if (fNormalCode == 0)	// normal is [1,0,0]
	{
		temp_bonds2 = temp_bonds;
		fBonds = temp_bonds2;
		fBonds *= -1.0;
	}
	else if (fNormalCode == 1)
		fBonds = temp_bonds;	// this table is the default, i.e. [-1,0,0]
	else if (fNormalCode == 2)	// rotate [-1,0,0] to [0,1,0]
	{
		temp_bonds2 = temp_bonds;
		blah1 = RotationMatrixA(piby2);
		for (int i = 0; i < kEAMFCC3DSurfBonds; i++)
		{
			temp_bonds2.RowCopy(i,asdf);	// take bond
			blah1.Multx(asdf,prod);		// rotate bond via rotation matrix
			temp_bonds2.SetRow(i,prod);	// place new bond back into temp_bonds
		}
		fBonds = temp_bonds2;
	}
	else if (fNormalCode == 3)	// rotate [-1,0,0] to [0,-1,0]
	{
		temp_bonds2 = temp_bonds;
		blah1 = RotationMatrixA(-piby2);
		for (int i = 0; i < kEAMFCC3DSurfBonds; i++)
		{
			temp_bonds2.RowCopy(i,asdf);	// take bond
			blah1.Multx(asdf,prod);		// rotate bond via rotation matrix
			temp_bonds2.SetRow(i,prod);	// place new bond back into temp_bonds
		}	
		fBonds = temp_bonds2;
	}
	else if (fNormalCode == 4)	// rotate [-1,0,0] to [0,0,1]
	{
		temp_bonds2 = temp_bonds;
		blah1 = RotationMatrixB(-piby2);
		for (int i = 0; i < kEAMFCC3DSurfBonds; i++)
		{
			temp_bonds2.RowCopy(i,asdf);	// take bond
			blah1.Multx(asdf,prod);		// rotate bond via rotation matrix
			temp_bonds2.SetRow(i,prod);	// place new bond back into temp_bonds
		}	
		fBonds = temp_bonds2;
	}	
	else if (fNormalCode == 5)	// rotate [-1,0,0] to [0,0,-1]
	{
		temp_bonds2 = temp_bonds;
		blah1 = RotationMatrixB(piby2);
		for (int i = 0; i < kEAMFCC3DSurfBonds; i++)
		{
			temp_bonds2.RowCopy(i,asdf);	// take bond
			blah1.Multx(asdf,prod);		// rotate bond via rotation matrix
			temp_bonds2.SetRow(i,prod);	// place new bond back into temp_bonds
		}	
		fBonds = temp_bonds2;
	}	
	
	/* copy into reference table */
	/* AVOID HARD CODING NUMBER OF BONDS DUE TO {100} SURFACES */
	/* OBVIATED BY fBonds ABOVE
	for (int i = 0; i < 78; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			fBonds(i,j) = bonddata[i][j];
	
	*/
	/* scale to correct lattice parameter */				     		
	fBonds *= fLatticeParameter;
	fBulkBonds *= fLatticeParameter;
	fSurf1Bonds *= fLatticeParameter;
	fSurf2Bonds *= fLatticeParameter;
}

/*************************************************************************
 * Private
 *************************************************************************/
 
 /* Rotate bonds with [-1,0,0] normal to bonds with [0,1,0]-type normals */
dMatrixT EAMFCC3DSym_surf::RotationMatrixA(const double angle)
 {
	dMatrixT rmatrix(3);
	rmatrix = 0.0;
    rmatrix(0,0) = cos(angle);
	rmatrix(0,1) = sin(angle);
	rmatrix(1,0) = -sin(angle);
	rmatrix(1,1) = cos(angle);
	rmatrix(2,2) = 1.0;
	
	return rmatrix;
 }
 
/* Rotate bonds with [-1,0,0] normal to bonds with [0,0,1]-type normals */
dMatrixT EAMFCC3DSym_surf::RotationMatrixB(const double angle)
{
	dMatrixT rmatrix(3);
	rmatrix = 0.0;
    rmatrix(0,0) = cos(angle);
	rmatrix(0,2) = -sin(angle);
	rmatrix(1,1) = 1.0;
	rmatrix(2,0) = sin(angle);
	rmatrix(2,2) = cos(angle);
	
	return rmatrix;
}
