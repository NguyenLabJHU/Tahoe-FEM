/* $Id: EAMFCC3DSym_surf.cpp,v 1.3 2006-06-03 23:03:45 hspark Exp $ */
/* created: paklein (12/06/1996) */
#include "EAMFCC3DSym_surf.h"

/* bond parameters */
const int kEAMFCC3DSymNumBonds = 78;	// changed from 27 for symmetric & bulk
const int kEAMFCC3DNumLatticeDim =  3;

using namespace Tahoe;

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
	/* dimension work space */
	fBondCounts.Dimension(kEAMFCC3DSymNumBonds);
	fDefLength.Dimension(kEAMFCC3DSymNumBonds);
	fBonds.Dimension(kEAMFCC3DSymNumBonds, kEAMFCC3DNumLatticeDim);

	dArray2DT temp_bonds, temp_bonds2;
	temp_bonds.Dimension(kEAMFCC3DSymNumBonds, 3);	// temporary bond table before rotation
	temp_bonds2.Dimension(kEAMFCC3DSymNumBonds, 3);	// Currently have # of bonds for {100} surfaces

	/* all bonds appear ONCE (change from regular EAMFCC3DSym) */
	fBondCounts = 1;
	
	/* clear deformed lengths for now */
	fDefLength = 0.0;

	/* New bond table for surface clusters - change dimensions! */
	/* AVOID HARD CODING NUMBER OF BONDS SPECIFIC TO {100} */
	double bonddata[kEAMFCC3DSymNumBonds][kEAMFCC3DNumLatticeDim] = {
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
	
	/* Rotate Bond Tables based on fNormalCode and rotation matrices */
	/* Create temporary bond table temp_bonds that combines bonddata */
	for (int i = 0; i < kEAMFCC3DSymNumBonds; i++)
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
		for (int i = 0; i < kEAMFCC3DSymNumBonds; i++)
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
		for (int i = 0; i < kEAMFCC3DSymNumBonds; i++)
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
		for (int i = 0; i < kEAMFCC3DSymNumBonds; i++)
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
		for (int i = 0; i < kEAMFCC3DSymNumBonds; i++)
		{
			temp_bonds2.RowCopy(i,asdf);	// take bond
			blah1.Multx(asdf,prod);		// rotate bond via rotation matrix
			temp_bonds2.SetRow(i,prod);	// place new bond back into temp_bonds
		}	
		fBonds = temp_bonds2;
	}
	
	/* copy into reference table */
	/* OBVIATED BY fBonds ABOVE */
	//for (int i = 0; i < kEAMFCC3DSymNumBonds; i++)
	//	for (int j = 0; j < 3; j++)
	//		fBonds(i,j) = bonddata[i][j];
			
	/* scale to correct lattice parameter */				     		
	fBonds *= fLatticeParameter;
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
