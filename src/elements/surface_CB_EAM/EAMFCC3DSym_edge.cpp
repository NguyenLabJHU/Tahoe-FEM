/* $Id: EAMFCC3DSym_edge.cpp,v 1.3 2009-06-04 04:28:00 hspark Exp $ */
/* created: paklein (12/06/1996) */
#include "EAMFCC3DSym_edge.h"

using namespace Tahoe;

/* Bond table parameters */
const int kEAMFCC3DEdgeBonds        = 412;	// total sum of neighbors
const int kEAMFCC3DEdgeBonds1        = 15;	// 15 atoms for type-1 (edge) unit cell x 2 atoms = 30
const int kEAMFCC3DEdgeBonds2        = 22;	// 22 atoms in type-2 unit cell x 2 atoms = 44
const int kEAMFCC3DEdgeBonds3        = 25;	// 25 atoms in type-3 unit cell x 1 atom = 25
const int kEAMFCC3DEdgeBonds4        = 22;	// 22 atoms in type-4 unit cell x 2 atoms = 44
const int kEAMFCC3DEdgeBonds5        = 25;	// 25 atoms in type-5 unit cell x 1 atom = 25
const int kEAMFCC3DEdgeBonds6        = 32;	// 32 atoms in type-6 unit cell x 3 atoms = 96
const int kEAMFCC3DEdgeBonds7        = 37;	// 37 atoms in type-7 unit cell x 2 atoms = 74
const int kEAMFCC3DEdgeBonds8        = 37;	// 37 atoms in type-8 unit cell x 2 atoms = 74
const int kEAMFCC3DNumLatticeDim 	=  3;
const int kEAMFCC3DNumAtomsPerCell	=  4;
const int kEAMFCC3DNumAtomsPerArea  =  2;
const double piby2 = 4.0 * atan(1.0) / 2.0;

/* constructor */
EAMFCC3DSym_edge::EAMFCC3DSym_edge(int nshells, int normal):
	EAMFCC3D_edge(nshells, normal),
	fNormalCode(normal)
{
	SetName("FCC_EAM_Cauchy-Born");
}

/**********************************************************************
 * Protected
 **********************************************************************/
	
void EAMFCC3DSym_edge::LoadBondTable(void)
{
	/* dimension work space - ARE THESE DIMENSIONS CORRECT? */
	/* RENAME VARIABLES IN BondTableT!!! */
	fBondCounts.Dimension(kEAMFCC3DEdgeBonds);
	fEdgeCounts.Dimension(kEAMFCC3DEdgeBonds);
	fDefLength.Dimension(kEAMFCC3DEdgeBonds);
	fDefEdge.Dimension(kEAMFCC3DEdgeBonds);
	fBonds.Dimension(kEAMFCC3DEdgeBonds, kEAMFCC3DNumLatticeDim);
	fEdgeBonds.Dimension(kEAMFCC3DEdgeBonds,3);
	fEdgeType.Dimension(kEAMFCC3DEdgeBonds);

	dArray2DT temp_edge, temp_edge1, temp_edge2, tempedge;
	temp_edge.Dimension(kEAMFCC3DEdgeBonds, 3);
	temp_edge1.Dimension(kEAMFCC3DEdgeBonds, 3);
	temp_edge2.Dimension(kEAMFCC3DEdgeBonds, 3);
	tempedge.Dimension(kEAMFCC3DEdgeBonds, 3);

	/* all bonds appear once */
	fBondCounts = 1;
	fEdgeCounts = 1;
	
	/* clear deformed lengths for now */
	fDefLength = 0.0;
	fDefEdge = 0.0;

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

	/* Bond table for an atom on the edge */
	double edgebond[kEAMFCC3DEdgeBonds1][kEAMFCC3DNumLatticeDim] = {
		{0.5, 0.5, 0.0}, // Surface cluster (8 nearest neighbors)
		{0.5, -0.5, 0.0},
		{0.5, 0.0, -0.5},
		{0.0, 0.5, -0.5},
		{0.0, -0.5, -0.5},
		{1.0, 0.0, 0.0}, // Surface cluster (5 2nd shell neighbors)
		{0.0, 1.0, 0.0},
		{0.0, -1.0, 0.0},
		{0.0, 0.0, -1.0},
		{1.0, 0.5, -0.5},
		{0.5, 1.0, -0.5},
		{0.5, 0.5, -1.0},
		{1.0, -0.5, -0.5},
		{0.5, -1.0, -0.5},
		{0.5, -0.5, -1.0}
	};
	/* temporary edgebond interaction map */
	// edgebondmap = {1,1,3,1,1,2,0,0,2,4,3,4,4,3,4};

	/* Bond table for z=-0.5, x=0, y=variable (2 atoms) */
	double edgebond2[kEAMFCC3DEdgeBonds2][kEAMFCC3DNumLatticeDim] = {
		{0, 0, -1.}, // 5
		{0, -1., 0}, // 6
		{0, 1., 0}, // 6
		{1., 0, 0}, // 8
		{0.5, 0, -0.5}, // 8
		{0.5, 0, 0.5}, // 6
		{0, -0.5, -0.5}, // 5
		{0, -0.5, 0.5}, // 1
		{0, 0.5, -0.5}, // 5
		{0, 0.5, 0.5}, // 1
		{0.5, -0.5, 0}, // 7
		{0.5, 0.5, 0}, // 7
		{0.5, -0.5, -1.}, // 8
		{0.5, 0.5, -1.}, // 8
		{0.5, -1., -0.5}, // 8
		{0.5, -1., 0.5}, // 6
		{0.5, 1., -0.5}, // 8
		{0.5, 1., 0.5}, // 6
		{1., -0.5, -0.5}, // 9
		{1., -0.5, 0.5}, // 5
		{1., 0.5, -0.5}, // 9
		{1., 0.5, 0.5} // 5
	};	
	/* temporary edgebond1 interaction map */
	// edgebond1map = {5,6,6,8,8,6,5,1,5,1,7,7,8,8,8,6,8,6,9,5,9,5};

	/* Bond table for z=-1.0, x=0, y=variable (1 atom) */
	double edgebond3[kEAMFCC3DEdgeBonds3][kEAMFCC3DNumLatticeDim] = {
		{0, 0, -1.}, // 10
		{0, 0, 1.}, // 2
		{0, -1., 0}, // 10
		{0, 1., 0}, // 10
		{1., 0, 0}, // 13
		{0.5, 0, -0.5}, // 12
		{0.5, 0, 0.5}, // 11
		{0, -0.5, -0.5}, // 10
		{0, -0.5, 0.5}, // 5
		{0, 0.5, -0.5}, // 10
		{0, 0.5, 0.5}, // 5
		{0.5, -0.5, 0}, // 12
		{0.5, 0.5, 0}, // 12
		{0.5, -0.5, -1.}, // 12
		{0.5, -0.5, 1.}, // 5 
		{0.5, 0.5, -1.}, // 12
		{0.5, 0.5, 1.}, // 5
		{0.5, -1., -0.5}, // 12
		{0.5, -1., 0.5}, // 11
		{0.5, 1., -0.5}, // 12
		{0.5, 1., 0.5}, // 11
		{1., -0.5, -0.5}, // 13
		{1., -0.5, 0.5}, // 12
		{1., 0.5, -0.5}, // 13
		{1., 0.5, 0.5} // 12
	};	
	/* temporary edgebond3 interaction map */
	// edgebond2map = {10,2,10,10,13,12,11,10,5,10,5,12,12,12,5,12,5,12,11,12,11,13,12,13,12}
	
	/* Bond table for z=0, x=0.5, y=variable (2 atoms) */
	double edgebond4[kEAMFCC3DEdgeBonds4][kEAMFCC3DNumLatticeDim] = {
		{0, 0, -1.}, // 8
		{0, -1., 0}, // 6
		{0, 1., 0}, // 6
		{1., 0, 0}, // 5
		{-0.5, 0, -0.5}, // 6
		{0.5, 0, -0.5}, // 8
		{0, -0.5, -0.5}, // 7
		{0, 0.5, -0.5}, // 7
		{-0.5, -0.5, 0}, // 1
		{-0.5, 0.5, 0}, // 1
		{0.5, -0.5, 0}, // 5
		{0.5, 0.5, 0}, // 5
		{-0.5, -0.5, -1.}, // 5
		{-0.5, 0.5, -1.}, // 5
		{0.5, -0.5, -1.}, // 9
		{0.5, 0.5, -1.}, // 9
		{-0.5, -1., -0.5}, // 6
		{-0.5, 1., -0.5}, // 6
		{0.5, -1., -0.5}, // 8
		{0.5, 1., -0.5}, // 8
		{1., -0.5, -0.5}, // 8
		{1., 0.5, -0.5} // 8
	};	
	/* temporary edgebond4 interaction map */
	// edgebond4map = {8,6,6,5,6,8,7,7,1,1,5,5,5,5,9,9,6,6,8,8,8,8};

	/* Bond table for z=0, x=1.0, y=variable (1 atom) */
	double edgebond5[kEAMFCC3DEdgeBonds5][kEAMFCC3DNumLatticeDim] = {
		{0, 0, -1.}, // 13
		{0, -1., 0}, // 10
		{0, 1., 0}, // 10
		{-1., 0, 0}, // 2
		{1., 0, 0}, // 10
		{-0.5, 0, -0.5}, // 11
		{0.5, 0, -0.5}, // 12
		{0, -0.5, -0.5}, // 12
		{0, 0.5, -0.5}, // 12
		{-0.5, -0.5, 0}, // 5
		{-0.5, 0.5, 0}, // 5
		{0.5, -0.5, 0}, // 10
		{0.5, 0.5, 0}, // 10
		{-0.5, -0.5, -1.}, // 12
		{-0.5, 0.5, -1.}, // 12
		{0.5, -0.5, -1.}, // 13
		{0.5, 0.5, -1.}, // 13
		{-0.5, -1., -0.5}, // 11 
		{-0.5, 1., -0.5}, // 11
		{0.5, -1., -0.5}, // 12
		{0.5, 1., -0.5}, // 12
		{-1., -0.5, -0.5}, // 5
		{-1., 0.5, -0.5}, // 5
		{1., -0.5, -0.5}, // 12
		{1., 0.5, -0.5} // 12
	};	
	/* temporary edgebond5 interaction map */
	// edgebond5map = {13,10,10,2,10,11,12,12,12,5,5,10,10,12,12,13,13,11,11,12,12,5,5,12,12};

	/* Bond table for x=0.5, z=-0.5, y=variable (3 atoms) (quasi surface 2) */
	double edgebond6[kEAMFCC3DEdgeBonds6][kEAMFCC3DNumLatticeDim] = {
		{0, 0, -1.}, // 15 
		{0, -1., 0}, // 14
		{0, 1., 0}, // 14
		{1., 0, 0}, // 15
		{-0.5, 0, -0.5}, // 11
		{-0.5, 0, 0.5}, // 3
		{0.5, 0, -0.5}, // 16
		{0.5, 0, 0.5}, // 11
		{0, -0.5, -0.5}, // 15
		{0, -0.5, 0.5}, // 7
		{0, 0.5, -0.5}, // 15
		{0, 0.5, 0.5}, // 7
		{-0.5, -0.5, 0}, // 3
		{-0.5, 0.5, 0}, // 3
		{0.5, -0.5, 0}, // 15
		{0.5, 0.5, 0}, // 15
		{-0.5, -0.5, -1.}, // 11
		{-0.5, 0.5, -1.}, // 11
		{0.5, -0.5, -1.}, // 16
		{0.5, 0.5, -1.}, // 16
		{-0.5, -1., -0.5}, // 11
		{-0.5, -1., 0.5}, // 3
		{-0.5, 1., -0.5}, // 11
		{-0.5, 1., 0.5}, // 3
		{0.5, -1., -0.5}, // 16
		{0.5, -1., 0.5}, // 11
		{0.5, 1., -0.5}, // 16
		{0.5, 1., 0.5}, // 11
		{1., -0.5, -0.5}, // 16
		{1., -0.5, 0.5}, // 11
		{1., 0.5, -0.5}, // 16
		{1., 0.5, 0.5} // 11
	};	
	/* temporary edgebond6 interaction map */
	// edgebond6map = {15,14,14,15,11,3,16,11,15,7,15,7,3,3,15,15,11,11,16,16,11,3,11,3,16,11,16,11,16,11,16,11}

	/* Bond table for x=0.5, z=-1.0, y=variable (2 atoms) */
	double edgebond7[kEAMFCC3DEdgeBonds7][kEAMFCC3DNumLatticeDim] = {
		{0, 0, -1.}, // 18
		{0, 0, 1.}, // 8
		{0, -1., 0}, // 18
		{0, 1., 0}, // 18
		{1., 0, 0}, // 17
		{-0.5, 0, -0.5}, // 12
		{-0.5, 0, 0.5}, // 8
		{0.5, 0, -0.5}, // 17
		{0.5, 0, 0.5}, // 18
		{0, -0.5, -0.5}, // 18
		{0, -0.5, 0.5}, // 15
		{0, 0.5, -0.5}, // 18
		{0, 0.5, 0.5}, // 15
		{-0.5, -0.5, 0}, // 12
		{-0.5, 0.5, 0}, // 12
		{0.5, -0.5, 0}, // 17
		{0.5, 0.5, 0}, // 17
		{-0.5, -0.5, -1.}, // 12
		{-0.5, -0.5, 1.}, // 4
		{-0.5, 0.5, -1.}, // 12
		{-0.5, 0.5, 1.}, // 4
		{0.5, -0.5, -1.}, // 17
		{0.5, -0.5, 1.}, // 12
		{0.5, 0.5, -1.}, // 17
		{0.5, 0.5, 1.}, // 12
		{-0.5, -1., -0.5}, // 12
		{-0.5, -1., 0.5}, // 8
		{-0.5, 1., -0.5}, // 12
		{-0.5, 1., 0.5}, // 8
		{0.5, -1., -0.5}, // 17
		{0.5, -1., 0.5}, // 18
		{0.5, 1., -0.5}, // 17
		{0.5, 1., 0.5}, // 18
		{1., -0.5, -0.5}, // 17
		{1., -0.5, 0.5}, // 18
		{1., 0.5, -0.5}, // 17
		{1., 0.5, 0.5} // 18
	};	
	/* temporary edgebond7 interaction map */
	// edgebond7map = {18,8,18,18,17,12,8,17,18,18,15,18,15,12,12,17,17,12,4,12,4,17,12,17,12,12,8,12,8,17,18,17,18,17,18,17,18}

	/* Bond table for x=1.0, z=-0.5, y=variable (2 atoms) */
	double edgebond8[kEAMFCC3DEdgeBonds8][kEAMFCC3DNumLatticeDim] = {
		{0, 0, -1.}, // 17
		{0, -1., 0}, // 18
		{0, 1., 0}, // 18
		{-1., 0, 0}, // 8
		{1., 0, 0}, // 18
		{-0.5, 0, -0.5}, // 18
		{-0.5, 0, 0.5}, // 8
		{0.5, 0, -0.5}, // 17
		{0.5, 0, 0.5}, // 12
		{0, -0.5, -0.5}, // 17
		{0, -0.5, 0.5}, // 12
		{0, 0.5, -0.5}, // 17
		{0, 0.5, 0.5}, // 12
		{-0.5, -0.5, 0}, // 15
		{-0.5, 0.5, 0}, // 15
		{0.5, -0.5, 0}, // 18
		{0.5, 0.5, 0}, // 18
		{-0.5, -0.5, -1.}, // 18
		{-0.5, 0.5, -1.}, // 18
		{0.5, -0.5, -1.}, // 17
		{0.5, 0.5, -1.}, // 17
		{-0.5, -1., -0.5}, // 18
		{-0.5, -1., 0.5}, // 8
		{-0.5, 1., -0.5}, // 18
		{-0.5, 1., 0.5}, // 8
		{0.5, -1., -0.5}, // 17
		{0.5, -1., 0.5}, // 12
		{0.5, 1., -0.5}, // 17
		{0.5, 1., 0.5}, // 12
		{-1., -0.5, -0.5}, // 12
		{-1., -0.5, 0.5}, // 4
		{-1., 0.5, -0.5}, // 12
		{-1., 0.5, 0.5}, // 4
		{1., -0.5, -0.5}, // 17
		{1., -0.5, 0.5}, // 12
		{1., 0.5, -0.5}, // 17
		{1., 0.5, 0.5} // 12
	};
	/* temporary edgebond8 interaction map */
	// edgebond8map = {17,18,18,8,18,18,8,17,12,17,12,17,12,15,15,18,18,18,18,17,17,18,8,18,8,17,12,17,12,12,4,12,4,17,12,17,12};

	/* Copy bond table into array */
	for (int i = 0; i < kEAMFCC3DSurf1Bonds; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			temp_edge(i,j) = edgebond[i][j];

	/* IMPLEMENT INTERACTION TABLE HERE FOR EDGE ATOM FOR CORRECT ELECTRON DENSITIES FOR
		OTHER EDGE/SURFACE ATOMS */
	/* work space arrays for storing interaction types */
	iArrayT allbonds(412);
	/* Interaction type key */
	// 0 = edge/edge
	// 1 = edge/1 atom away on surface 1
	// 2 = edge/2 atoms away on surface 1
	// 3 = edge/quasi surface 2 (32 neighbors)
	// 4 = edge/real surface 2 (37 neighbors)
	// 5 = 1 atom away on surface 1 / 2 atoms away on surface 1
	// 6 = 1 atom away on surface 1 / 1 atom away on surface 1
	// 7 = 1 atom away on surface 1 / quasi surface 2
	// 8 = 1 atom away on surface 1 / real surface 2
	// 9 = 1 atom away on surface 1 / bulk
	//10 = 2 atoms away on surface 1 / 2 atoms away on surface 1
	//11 = 2 atoms away on surface 1 / quasi surface 2
	//12 = 2 atoms away on surface 1 / real surface 2
	//13 = 2 atoms away on surface 1 / bulk
	//14 = quasi surface 2 / quasi surface 2
	//15 = quasi surface 2 / real surface 2
	//16 = quasi surface 2 / bulk
	//17 = real surface 2 / bulk
	//18 = real surface 2 / real surface 2
//	int edge1n[412]={};
//	for (int i = 0; i < kEAMFCC3DEdgeBonds; i++)
//		allbonds[i] = edge1n[i];
	
	//fEdgeType.CopyIn(0, allbonds);
	
	/* New bond table for surface clusters -  */
	double bonddata[kEAMFCC3DEdgeBonds][kEAMFCC3DNumLatticeDim] = {
		{0.5, 0.5, 0.0}, // Surface cluster (8 nearest neighbors)
		{0.5, -0.5, 0.0},
		{0.5, 0.0, -0.5},
		{0.0, 0.5, -0.5},
		{0.0, -0.5, -0.5},
		{1.0, 0.0, 0.0}, // Surface cluster (5 2nd shell neighbors)
		{0.0, 1.0, 0.0},
		{0.0, -1.0, 0.0},
		{0.0, 0.0, -1.0},
		{1.0, 0.5, -0.5},
		{0.5, 1.0, -0.5},
		{0.5, 0.5, -1.0},
		{1.0, -0.5, -0.5},
		{0.5, -1.0, -0.5},
		{0.5, -0.5, -1.0}
	};
	
	/* Rotate Bond Tables based on fNormalCode and rotation matrices */
	/* Create temporary bond table temp_bonds that combines bonddata */
	for (int i = 0; i < kEAMFCC3DEdgeBonds; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			temp_edge1(i,j) = bonddata[i][j];
	
	/* Now manipulate temp_bonds */
	dMatrixT blah1(3);
	dArrayT asdf(3), prod(3), asdfb(3), prodb(3), asdfs1(3), prods1(3), asdfs2(3), prods2(3);
	if (fNormalCode == 0)	// normal is [1,0,0]
	{
		temp_edge2 = temp_edge1;
		fBonds = temp_edge2;
		fBonds *= -1.0;
		fEdgeBonds = temp_edge2;
		fEdgeBonds *= -1.0;
	}
	else if (fNormalCode == 1)
	{
		fBonds = temp_edge1;	// this table is the default, i.e. [-1,0,0]
		fEdgeBonds = temp_edge;
	}
	else if (fNormalCode == 2)	// rotate [-1,0,0] to [0,1,0]
	{
		temp_edge2 = temp_edge1;
		tempedge = temp_edge;
		blah1 = RotationMatrixA(piby2);
		for (int i = 0; i < kEAMFCC3DEdgeBonds; i++)
		{
			temp_edge2.RowCopy(i,asdf);	// take bond
			blah1.Multx(asdf,prod);		// rotate bond via rotation matrix
			temp_edge2.SetRow(i,prod);	// place new bond back into temp_bonds
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds; i++)
		{
			tempedge.RowCopy(i,asdfs1);
			blah1.Multx(asdfs1,prods1);
			tempedge.SetRow(i,prods1);
		}
		fBonds = temp_edge2;
		fEdgeBonds = tempedge;
	}
	else if (fNormalCode == 3)	// rotate [-1,0,0] to [0,-1,0]
	{
		temp_edge2 = temp_edge1;
		tempedge = temp_edge;
		blah1 = RotationMatrixA(-piby2);
		for (int i = 0; i < kEAMFCC3DEdgeBonds; i++)
		{
			temp_edge2.RowCopy(i,asdf);	// take bond
			blah1.Multx(asdf,prod);		// rotate bond via rotation matrix
			temp_edge2.SetRow(i,prod);	// place new bond back into temp_bonds
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds; i++)
		{
			tempedge.RowCopy(i,asdfs1);
			blah1.Multx(asdfs1,prods1);
			tempedge.SetRow(i,prods1);
		}
		fBonds = temp_edge2;
		fEdgeBonds = tempedge;
	}
	else if (fNormalCode == 4)	// rotate [-1,0,0] to [0,0,1]
	{
		temp_edge2 = temp_edge1;
		tempedge = temp_edge;
		blah1 = RotationMatrixB(-piby2);
		for (int i = 0; i < kEAMFCC3DEdgeBonds; i++)
		{
			temp_edge2.RowCopy(i,asdf);	// take bond
			blah1.Multx(asdf,prod);		// rotate bond via rotation matrix
			temp_edge2.SetRow(i,prod);	// place new bond back into temp_bonds
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds; i++)
		{
			tempedge.RowCopy(i,asdfs1);
			blah1.Multx(asdfs1,prods1);
			tempedge.SetRow(i,prods1);
		}
		fBonds = temp_edge2;
		fEdgeBonds = tempedge;
	}	
	else if (fNormalCode == 5)	// rotate [-1,0,0] to [0,0,-1]
	{
		temp_edge2 = temp_edge1;
		tempedge = temp_edge;
		blah1 = RotationMatrixB(piby2);
		for (int i = 0; i < kEAMFCC3DEdgeBonds; i++)
		{
			temp_edge2.RowCopy(i,asdf);	// take bond
			blah1.Multx(asdf,prod);		// rotate bond via rotation matrix
			temp_edge2.SetRow(i,prod);	// place new bond back into temp_bonds
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds; i++)
		{
			tempedge.RowCopy(i,asdfs1);
			blah1.Multx(asdfs1,prods1);
			tempedge.SetRow(i,prods1);
		}
		fBonds = temp_edge2;
		fEdgeBonds = tempedge;
	}	

	/* scale to correct lattice parameter */				     		
	fBonds *= fLatticeParameter;
	fEdgeBonds *= fLatticeParameter;
}

/*************************************************************************
 * Private
 *************************************************************************/
 
 /* Rotate bonds with [-1,0,0] normal to bonds with [0,1,0]-type normals */
dMatrixT EAMFCC3DSym_edge::RotationMatrixA(const double angle)
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
dMatrixT EAMFCC3DSym_edge::RotationMatrixB(const double angle)
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
