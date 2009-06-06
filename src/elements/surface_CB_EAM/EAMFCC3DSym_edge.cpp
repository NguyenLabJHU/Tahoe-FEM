/* $Id: EAMFCC3DSym_edge.cpp,v 1.7 2009-06-06 17:28:31 hspark Exp $ */
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
const int kEAMFCC3DBulkBonds         = 42;	// 42 atoms in a bulk unit cell
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
	/* RENAME VARIABLES IN BondTableT!!! */
//	fEdgeBondCounts.Dimension(kEAMFCC3DEdgeBonds); // all edge bonds (412)
	fBondCounts.Dimension(kEAMFCC3DEdgeBonds);	// all edge bonds (412)
 	fBulkBondCounts.Dimension(kEAMFCC3DBulkBonds); // all bulk bonds (42)
 	fEdge1Counts.Dimension(kEAMFCC3DEdgeBonds1);
	fEdge2Counts.Dimension(kEAMFCC3DEdgeBonds2);
	fEdge3Counts.Dimension(kEAMFCC3DEdgeBonds3);
	fEdge4Counts.Dimension(kEAMFCC3DEdgeBonds4);
	fEdge5Counts.Dimension(kEAMFCC3DEdgeBonds5);
	fEdge6Counts.Dimension(kEAMFCC3DEdgeBonds6);
	fEdge7Counts.Dimension(kEAMFCC3DEdgeBonds7);
	fEdge8Counts.Dimension(kEAMFCC3DEdgeBonds8);	
//	fDefEdgeLength.Dimension(kEAMFCC3DEdgeBonds);	// all deformed edge lengths (412)
	fDefLength.Dimension(kEAMFCC3DEdgeBonds);
	fDefBulkLength.Dimension(kEAMFCC3DBulkBonds);
 	fDefEdge1.Dimension(kEAMFCC3DEdgeBonds1);
	fDefEdge2.Dimension(kEAMFCC3DEdgeBonds2);
	fDefEdge3.Dimension(kEAMFCC3DEdgeBonds3);
	fDefEdge4.Dimension(kEAMFCC3DEdgeBonds4);
	fDefEdge5.Dimension(kEAMFCC3DEdgeBonds5);
	fDefEdge6.Dimension(kEAMFCC3DEdgeBonds6);
	fDefEdge7.Dimension(kEAMFCC3DEdgeBonds7);
	fDefEdge8.Dimension(kEAMFCC3DEdgeBonds8);	
//	fAllEdgeBonds.Dimension(kEAMFCC3DEdgeBonds,3);
	fBonds.Dimension(kEAMFCC3DEdgeBonds,3);
	fAllBulkBonds.Dimension(kEAMFCC3DBulkBonds,3);
 	fEdgeBonds1.Dimension(kEAMFCC3DEdgeBonds1,3);
	fEdgeBonds2.Dimension(kEAMFCC3DEdgeBonds2,3);
	fEdgeBonds3.Dimension(kEAMFCC3DEdgeBonds3,3);
	fEdgeBonds4.Dimension(kEAMFCC3DEdgeBonds4,3);
	fEdgeBonds5.Dimension(kEAMFCC3DEdgeBonds5,3);
	fEdgeBonds6.Dimension(kEAMFCC3DEdgeBonds6,3);
	fEdgeBonds7.Dimension(kEAMFCC3DEdgeBonds7,3);
	fEdgeBonds8.Dimension(kEAMFCC3DEdgeBonds8,3);
	fEdgeType.Dimension(kEAMFCC3DEdgeBonds);

	dArray2DT temp_edge1, temp_edge2, temp_edge3, temp_edge4, temp_edge5, temp_edge6, temp_edge7, temp_edge8, temp_bulk; 
	dArray2DT tempedge;
	temp_edge1.Dimension(kEAMFCC3DEdgeBonds1, 3);
	temp_edge2.Dimension(kEAMFCC3DEdgeBonds2, 3);
	temp_edge3.Dimension(kEAMFCC3DEdgeBonds3, 3);
	temp_edge4.Dimension(kEAMFCC3DEdgeBonds4, 3);
	temp_edge5.Dimension(kEAMFCC3DEdgeBonds5, 3);
	temp_edge6.Dimension(kEAMFCC3DEdgeBonds6, 3);
	temp_edge7.Dimension(kEAMFCC3DEdgeBonds7, 3);
	temp_edge8.Dimension(kEAMFCC3DEdgeBonds8, 3);
	temp_bulk.Dimension(kEAMFCC3DBulkBonds, 3);
	tempedge.Dimension(kEAMFCC3DEdgeBonds, 3);

	/* all bonds appear once */
//	fEdgeBondCounts = 1;
	fBondCounts = 1;
	fBulkBondCounts = 1;
 	fEdge1Counts = 1;
	fEdge2Counts = 1;
	fEdge3Counts = 1;
	fEdge4Counts = 1;
	fEdge5Counts = 1;
	fEdge6Counts = 1;
	fEdge7Counts = 1;
	fEdge8Counts = 1;
	
	/* clear deformed lengths for now */
//	fDefEdgeLength = 0.0;
	fDefLength = 0.0;
	fDefBulkLength = 0.0;
 	fDefEdge1 = 0.0;
	fDefEdge2 = 0.0;
	fDefEdge3 = 0.0;
	fDefEdge4 = 0.0;
	fDefEdge5 = 0.0;
	fDefEdge6 = 0.0;
	fDefEdge7 = 0.0;
	fDefEdge8 = 0.0;	

	/* undeformed bond data for bulk atom with 4th neighbor interactions */
	double bulkbond[kEAMFCC3DBulkBonds][kEAMFCC3DNumLatticeDim] = {
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

	/* Copy bond table into array */
	for (int i = 0; i < kEAMFCC3DBulkBonds; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			temp_bulk(i,j) = bulkbond[i][j];

	/* Bond table for an atom on the edge (2 atoms) */
	double edgebond1[kEAMFCC3DEdgeBonds1][kEAMFCC3DNumLatticeDim] = {
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
	
	/* Copy bond table into array */
	for (int i = 0; i < kEAMFCC3DEdgeBonds1; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			temp_edge1(i,j) = edgebond1[i][j];
			
	/* temporary edgebond1 interaction map */
	// edgebond1map = {1,1,3,1,1,2,0,0,2,4,3,4,4,3,4}; - WRONG
	// edgebond1map = {1,1,2,3,3,4,0,0,5,6,2,7,6,2,7}; - CORRECTED

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
	
	/* Copy bond table into array */
	for (int i = 0; i < kEAMFCC3DEdgeBonds2; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			temp_edge2(i,j) = edgebond2[i][j];
			
	/* temporary edgebond2 interaction map */
	// edgebond2map = {5,6,6,8,8,6,5,1,5,1,7,7,8,8,8,6,8,6,9,5,9,5}; - WRONG
	// edgebond2map = {8,9,9,10,11,12,8,3,8,3,13,13,11,11,11,12,11,12,14,15,14,15} - CORRECTED

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
	
	/* Copy bond table into array */
	for (int i = 0; i < kEAMFCC3DEdgeBonds3; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			temp_edge3(i,j) = edgebond3[i][j];
			
	/* temporary edgebond3 interaction map */
	// edgebond3map = {10,2,10,10,13,12,11,10,5,10,5,12,12,12,5,12,5,12,11,12,11,13,12,13,12} - WRONG
	// edgebond3map = {16,5,16,5,17,18,19,16,20,16,20,18,18,18,21,18,21,18,19,18,19,17,22,17,22}; - CORRECT
	
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
	
	/* Copy bond table into array */
	for (int i = 0; i < kEAMFCC3DEdgeBonds4; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			temp_edge4(i,j) = edgebond4[i][j];
			
	/* temporary edgebond4 interaction map */
	// edgebond4map = {8,6,6,5,6,8,7,7,1,1,5,5,5,5,9,9,6,6,8,8,8,8}; - WRONG
	// edgebond4map = {23,24,24,25,12,26,27,27,1,1,25,25,28,28,29,29,12,12,26,26,26,26};

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
	
	/* Copy bond table into array */
	for (int i = 0; i < kEAMFCC3DEdgeBonds5; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			temp_edge5(i,j) = edgebond5[i][j];
			
	/* temporary edgebond5 interaction map */
	// edgebond5map = {13,10,10,2,10,11,12,12,12,5,5,10,10,12,12,13,13,11,11,12,12,5,5,12,12}; - WRONG
	// edgebond5map = {30,31,31,4,31,32,33,33,33,25,25,31,31,34,34,30,30,32,32,33,33,15,15,33,33};

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
	
	/* Copy bond table into array */
	for (int i = 0; i < kEAMFCC3DEdgeBonds6; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			temp_edge6(i,j) = edgebond6[i][j];
			
	/* temporary edgebond6 interaction map */
	// edgebond6map = {15,14,14,15,11,3,16,11,15,7,15,7,3,3,15,15,11,11,16,16,11,3,11,3,16,11,16,11,16,11,16,11} - WRONG
	// edgebond6map = {35,36,36,37,19,2,38,32,35,27,35,27,13,13,37,37,19,19,38,38,19,2,19,2,38,32,38,32,38,32,38,32};

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
	
	/* Copy bond table into array */
	for (int i = 0; i < kEAMFCC3DEdgeBonds7; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			temp_edge7(i,j) = edgebond7[i][j];
			
	/* temporary edgebond7 interaction map */
	// edgebond7map = {18,8,18,18,17,12,8,17,18,18,15,18,15,12,12,17,17,12,4,12,4,17,12,17,12,12,8,12,8,17,18,17,18,17,18,17,18} - WRONG
	// edgebond7map = {39,23,39,39,40,18,11,40,41,39,35,39,35,18,18,40,40,18,7,18,7,40,34,40,34,18,11,18,11,40,41,40,41,40,41,40,41}; - CORRECT

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
	
	/* Copy bond table into array */
	for (int i = 0; i < kEAMFCC3DEdgeBonds8; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			temp_edge8(i,j) = edgebond8[i][j];
			
	/* temporary edgebond8 interaction map */
	// edgebond8map = {17,18,18,8,18,18,8,17,12,17,12,17,12,15,15,18,18,18,18,17,17,18,8,18,8,17,12,17,12,12,4,12,4,17,12,17,12}; - WRONG
	// edgebond8map = {42,43,43,10,43,41,26,42,33,42,33,42,33,37,37,43,43,41,41,42,42,41,26,41,26,42,33,42,33,22,6,22,6,42,33,42,33};  CORRECT
	
	/* IMPLEMENT INTERACTION TABLE HERE FOR EDGE ATOM FOR CORRECT ELECTRON DENSITIES FOR
		OTHER EDGE/SURFACE ATOMS */
	/* Interaction type key */
// 0:  Edge-Edge
// 1:  Edge-1 atom away (+Z)
// 2:  Edge-quasi surface 2 (both ÐX and +Z)
// 3:  Edge-1 atom away (-X)
// 4:  Edge-2 atoms away (+Z)
// 5:  Edge-2 atoms away (-X)
// 6:  Edge-real surface 2 (+Z)
// 7:  Edge-real surface 2 (-X)
// 8:  1 atom away (-X) Ð 2 atoms away (-X)
// 9:  1 atom away (-X) Ð 1 atom away (-X)
// 10:  1 atom away (-X) Ð real surface 2 (+Z)
// 11:  1 atom away (-X) Ð real surface 2 (-X)
// 12:  1 atom away (-X) Ð 1 atom away (+Z)
// 13:  1 atom away (-X) Ð quasi surface 2 (both ÐX and +Z)
// 14:  1 atom away (-X) Ð bulk
// 15:  1 atom away (-X) Ð 2 atoms away (+Z)
// 16:  2 atoms away (-X) Ð 2 atoms away (-X)
// 17:  2 atoms away (-X) Ð bulk
// 18:  2 atoms away (-X) Ð real surface 2 (-X)
// 19:  2 atoms away (-X) Ð quasi surface 2 (both ÐX and +Z)
// 20:  2 atoms away (-X) Ð 1 atom away (-X)
// 21:  2 atoms away (-X) Ð 1 atom away (+Z)
// 22:  2 atoms away (-X) Ð real surface 2 (+Z)
// 23:  1 atom away (+Z) Ð real surface 2 (-X)
// 24:  1 atom away (+Z) Ð 1 atom away (+Z)
// 25:  1 atom away (+Z) Ð 2 atoms away (+Z)
// 26:  1 atom away (+Z) Ð real surface 2 (+Z)
// 27:  1 atom away (+Z) Ð quasi surface 2 (both ÐX and +Z)
// 28:  1 atom away (+Z) Ð 2 atoms away (-X)
// 29:  1 atom away (+Z) Ð bulk
// 30:  2 atoms away (+Z) Ð bulk
// 31:  2 atoms away (+Z) Ð 2 atoms away (+Z)
// 32:  2 atoms away (+Z) Ð quasi surface 2 (-X and +Z)
// 33:  2 atoms away (+Z) Ð real surface 2 (+Z)
// 34:  2 atoms away (+Z) Ð real surface 2 (-X)
// 35:  quasi surface 2 Ð real surface 2 (-X)
// 36:  quasi surface 2 Ð quasi surface 2
// 37:  quasi surface 2 Ð real surface 2 (+Z)
// 38:  quasi surface 2 Ð bulk
// 39:  real surface 2 (-X) Ð real surface 2 (-X)
// 40:  real surface 2 (-X) Ð bulk
// 41:  real surface 2 (-X) Ð real surface 2 (+Z)
// 42:  real surface 2 (+Z) Ð bulk
// 43:  real surface 2 (+Z) Ð real surface 2 (+Z)
	// edge0 = bulk
	// edge1 = edge
	// edge2 = 1 atom away on surface 1 (-X normal)
	// edge3 = 2 atoms away on surface 1 (-X normal)
	// edge4 = 1 atom away on surface 1 (+Z normal)
	// edge5 = 2 atoms away on surface 1 (+Z normal)
	// edge6 = quasi surface 2
	// edge7 = real surface 2 (-X normal)
	// edge8 = real surface 2 (+Z normal)
	int edge1n[412]={1,1,2,3,3,4,0,0,5,6,2,7,6,2,7,1,1,2,3,3,4,0,0,5,6,2,7,6,2,7,8,9,9,10,11,12,8,3,8,3,13,13,11,11,11,12,11,12,14,15,14,15,8,9,9,10,11,12,8,3,8,3,13,13,11,11,11,12,11,12,14,15,14,15,16,5,16,5,17,18,19,16,20,16,20,18,18,18,21,18,21,18,19,18,19,17,22,17,22,23,24,24,25,12,26,27,27,1,1,25,25,28,28,29,29,12,12,26,26,26,26,23,24,24,25,12,26,27,27,1,1,25,25,28,28,29,29,12,12,26,26,26,26,30,31,31,4,31,32,33,33,33,25,25,31,31,34,34,30,30,32,32,33,33,15,15,33,33,35,36,36,37,19,2,38,32,35,27,35,27,13,13,37,37,19,19,38,38,19,2,19,2,38,32,38,32,38,32,38,32,35,36,36,37,19,2,38,32,35,27,35,27,13,13,37,37,19,19,38,38,19,2,19,2,38,32,38,32,38,32,38,32,35,36,36,37,19,2,38,32,35,27,35,27,13,13,37,37,19,19,38,38,19,2,19,2,38,32,38,32,38,32,38,32,39,23,39,39,40,18,11,40,41,39,35,39,35,18,18,40,40,18,7,18,7,40,34,40,34,18,11,18,11,40,41,40,41,40,41,40,41,39,23,39,39,40,18,11,40,41,39,35,39,35,18,18,40,40,18,7,18,7,40,34,40,34,18,11,18,11,40,41,40,41,40,41,40,41,42,43,43,10,43,41,26,42,33,42,33,42,33,37,37,43,43,41,41,42,42,41,26,41,26,42,33,42,33,22,6,22,6,42,33,42,33,42,43,43,10,43,41,26,42,33,42,33,42,33,37,37,43,43,41,41,42,42,41,26,41,26,42,33,42,33,22,6,22,6,42,33,42,33};  		
	
	/* work space arrays for storing interaction types */
	iArrayT allbonds(412);
	for (int i = 0; i < kEAMFCC3DEdgeBonds; i++)
		allbonds[i] = edge1n[i];
	
	fEdgeType.CopyIn(0, allbonds);
	
	/* Bond table containing all 412 bonds */
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
		{0.5, -0.5, -1.0},
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
		{0.5, -0.5, -1.0}, // END OF 2 TYPE-1 ATOMS
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
		{1., 0.5, 0.5}, // 5		
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
		{1., 0.5, 0.5}, // 5 // END OF 2 TYPE-2 ATOMS		
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
		{1., 0.5, 0.5}, // 12 // END OF 1 TYPE-3 ATOM
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
		{1., 0.5, -0.5}, // 8
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
		{1., 0.5, -0.5}, // 8 // END OF 2 TYPE-4 ATOMS
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
		{1., 0.5, -0.5}, // 12 // END OF 1 TYPE-5 ATOM
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
		{1., 0.5, 0.5}, // 11		
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
		{1., 0.5, 0.5}, // 11
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
		{1., 0.5, 0.5}, // 11 // END OF 3 TYPE-6 ATOMS
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
		{1., 0.5, 0.5}, // 18
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
		{1., 0.5, 0.5}, // 18 // END OF 2 TYPE-7 ATOMS
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
		{1., 0.5, 0.5}, // 12
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
		{1., 0.5, 0.5} // 12 // END OF 2 TYPE-8 ATOMS	
	};
	
	/* Rotate Bond Tables based on fNormalCode and rotation matrices */
	/* Create temporary bond table temp_bonds that combines bonddata */
	for (int i = 0; i < kEAMFCC3DEdgeBonds; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			tempedge(i,j) = bonddata[i][j];
	
	/* Now manipulate temp_bonds */
	dMatrixT blah1(3);
	dArrayT asdfe(3), prode(3), asdfb(3), prodb(3), asdfe1(3), prode1(3), asdfe2(3), prode2(3);
	dArrayT asdfe3(3), prode3(3), asdfe4(3), prode4(3), asdfe5(3), prode5(3), asdfe6(3), prode6(3);
	dArrayT asdfe7(3), prode7(3), asdfe8(3), prode8(3);
	
	dArray2DT tedge2, tbulk, te1, te2, te3, te4, te5, te6, te7, te8;
	tedge2.Dimension(kEAMFCC3DEdgeBonds, 3);
	tbulk.Dimension(kEAMFCC3DBulkBonds, 3);
	te1.Dimension(kEAMFCC3DEdgeBonds1, 3);
	te2.Dimension(kEAMFCC3DEdgeBonds2, 3);
	te3.Dimension(kEAMFCC3DEdgeBonds3, 3);
	te4.Dimension(kEAMFCC3DEdgeBonds4, 3);
	te5.Dimension(kEAMFCC3DEdgeBonds5, 3);
	te6.Dimension(kEAMFCC3DEdgeBonds6, 3);
	te7.Dimension(kEAMFCC3DEdgeBonds7, 3);
	te8.Dimension(kEAMFCC3DEdgeBonds8, 3);
	
	if (fNormalCode == 0)	// normal is [1,0,0]
	{
		tedge2 = tempedge;
//		fAllEdgeBonds = tedge2;
//		fAllEdgeBonds *= -1.0;
		fBonds = tedge2;
		fBonds *= -1.0;
		fAllBulkBonds = temp_bulk;
		fAllBulkBonds *= -1.0;
 		fEdgeBonds1 = temp_edge1;
 		fEdgeBonds1 *= -1.0;
		fEdgeBonds2 = temp_edge2;
		fEdgeBonds2 *= -1.0;
		fEdgeBonds3 = temp_edge3;
		fEdgeBonds3 *= -1.0;
		fEdgeBonds4 = temp_edge4;
		fEdgeBonds4 *= -1.0;
		fEdgeBonds5 = temp_edge5;
		fEdgeBonds5 *= -1.0;
		fEdgeBonds6 = temp_edge6;
		fEdgeBonds6 *= -1.0;
		fEdgeBonds7 = temp_edge7;
		fEdgeBonds7 *= -1.0;
		fEdgeBonds8 = temp_edge8;
		fEdgeBonds8 *= -1.0;
	}
	else if (fNormalCode == 1) // this table is the default orientation, i.e. [-1,0,0]
	{
//		fAllEdgeBonds = tedge2;
		fBonds = tempedge;
		fAllBulkBonds = temp_bulk;
 		fEdgeBonds1 = temp_edge1;
		fEdgeBonds2 = temp_edge2;
		fEdgeBonds3 = temp_edge3;
		fEdgeBonds4 = temp_edge4;
		fEdgeBonds5 = temp_edge5;
		fEdgeBonds6 = temp_edge6;
		fEdgeBonds7 = temp_edge7;
		fEdgeBonds8 = temp_edge8;
	}
	else if (fNormalCode == 2)	// rotate [-1,0,0] to [0,1,0]
	{
		tedge2 = tempedge;
		tbulk = temp_bulk;
 		te1 = temp_edge1;
		te2 = temp_edge2;
		te3 = temp_edge3;
		te4 = temp_edge4;
		te5 = temp_edge5;
		te6 = temp_edge6;
		te7 = temp_edge7;
		te8 = temp_edge8;		
		blah1 = RotationMatrixA(piby2);
		for (int i = 0; i < kEAMFCC3DEdgeBonds; i++)
		{
			tedge2.RowCopy(i,asdfe);	// take bond
			blah1.Multx(asdfe,prode);		// rotate bond via rotation matrix
			tedge2.SetRow(i,prode);	// place new bond back into temp_bonds
		}
		for (int i = 0; i < kEAMFCC3DBulkBonds; i++)
		{
			tbulk.RowCopy(i,asdfb);
			blah1.Multx(asdfb,prodb);
			tbulk.SetRow(i,prodb);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds1; i++)
		{
			te1.RowCopy(i,asdfe1);
			blah1.Multx(asdfe1,prode1);
			te1.SetRow(i,prode1);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds2; i++)
		{
			te2.RowCopy(i,asdfe2);
			blah1.Multx(asdfe2,prode2);
			te2.SetRow(i,prode2);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds3; i++)
		{
			te3.RowCopy(i,asdfe3);
			blah1.Multx(asdfe3,prode3);
			te3.SetRow(i,prode3);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds4; i++)
		{
			te4.RowCopy(i,asdfe4);
			blah1.Multx(asdfe4,prode4);
			te4.SetRow(i,prode4);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds5; i++)
		{
			te5.RowCopy(i,asdfe5);
			blah1.Multx(asdfe5,prode5);
			te5.SetRow(i,prode5);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds6; i++)
		{
			te6.RowCopy(i,asdfe6);
			blah1.Multx(asdfe6,prode6);
			te6.SetRow(i,prode6);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds7; i++)
		{
			te7.RowCopy(i,asdfe7);
			blah1.Multx(asdfe7,prode7);
			te7.SetRow(i,prode7);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds8; i++)
		{
			te8.RowCopy(i,asdfe8);
			blah1.Multx(asdfe8,prode8);
			te8.SetRow(i,prode8);
		}
//		fAllEdgeBonds = tedge2;
		fBonds = tedge2;
		fAllBulkBonds = tbulk;
 		fEdgeBonds1 = te1;
		fEdgeBonds2 = te2;
		fEdgeBonds3 = te3;
		fEdgeBonds4 = te4;
		fEdgeBonds5 = te5;
		fEdgeBonds6 = te6;
		fEdgeBonds7 = te7;
		fEdgeBonds8 = te8;		
	}
	else if (fNormalCode == 3)	// rotate [-1,0,0] to [0,-1,0]
	{
		tedge2 = tempedge;
		tbulk = temp_bulk;
 		te1 = temp_edge1;
		te2 = temp_edge2;
		te3 = temp_edge3;
		te4 = temp_edge4;
		te5 = temp_edge5;
		te6 = temp_edge6;
		te7 = temp_edge7;
		te8 = temp_edge8;	
		blah1 = RotationMatrixA(-piby2);
		for (int i = 0; i < kEAMFCC3DEdgeBonds; i++)
		{
			tedge2.RowCopy(i,asdfe);	// take bond
			blah1.Multx(asdfe,prode);		// rotate bond via rotation matrix
			tedge2.SetRow(i,prode);	// place new bond back into temp_bonds
		}
		for (int i = 0; i < kEAMFCC3DBulkBonds; i++)
		{
			tbulk.RowCopy(i,asdfb);
			blah1.Multx(asdfb,prodb);
			tbulk.SetRow(i,prodb);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds1; i++)
		{
			te1.RowCopy(i,asdfe1);
			blah1.Multx(asdfe1,prode1);
			te1.SetRow(i,prode1);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds2; i++)
		{
			te2.RowCopy(i,asdfe2);
			blah1.Multx(asdfe2,prode2);
			te2.SetRow(i,prode2);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds3; i++)
		{
			te3.RowCopy(i,asdfe3);
			blah1.Multx(asdfe3,prode3);
			te3.SetRow(i,prode3);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds4; i++)
		{
			te4.RowCopy(i,asdfe4);
			blah1.Multx(asdfe4,prode4);
			te4.SetRow(i,prode4);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds5; i++)
		{
			te5.RowCopy(i,asdfe5);
			blah1.Multx(asdfe5,prode5);
			te5.SetRow(i,prode5);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds6; i++)
		{
			te6.RowCopy(i,asdfe6);
			blah1.Multx(asdfe6,prode6);
			te6.SetRow(i,prode6);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds7; i++)
		{
			te7.RowCopy(i,asdfe7);
			blah1.Multx(asdfe7,prode7);
			te7.SetRow(i,prode7);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds8; i++)
		{
			te8.RowCopy(i,asdfe8);
			blah1.Multx(asdfe8,prode8);
			te8.SetRow(i,prode8);
		}
//		fAllEdgeBonds = tedge2;
		fBonds = tedge2;
		fAllBulkBonds = tbulk;
		fEdgeBonds1 = te1;
		fEdgeBonds2 = te2;
		fEdgeBonds3 = te3;
		fEdgeBonds4 = te4;
		fEdgeBonds5 = te5;
		fEdgeBonds6 = te6;
		fEdgeBonds7 = te7;
		fEdgeBonds8 = te8;		
	}
	else if (fNormalCode == 4)	// rotate [-1,0,0] to [0,0,1]
	{
		tedge2 = tempedge;
		tbulk = temp_bulk;
 		te1 = temp_edge1;
		te2 = temp_edge2;
		te3 = temp_edge3;
		te4 = temp_edge4;
		te5 = temp_edge5;
		te6 = temp_edge6;
		te7 = temp_edge7;
		te8 = temp_edge8;	
		blah1 = RotationMatrixB(-piby2);
		for (int i = 0; i < kEAMFCC3DEdgeBonds; i++)
		{
			tedge2.RowCopy(i,asdfe);	// take bond
			blah1.Multx(asdfe,prode);		// rotate bond via rotation matrix
			tedge2.SetRow(i,prode);	// place new bond back into temp_bonds
		}
		for (int i = 0; i < kEAMFCC3DBulkBonds; i++)
		{
			tbulk.RowCopy(i,asdfb);
			blah1.Multx(asdfb,prodb);
			tbulk.SetRow(i,prodb);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds1; i++)
		{
			te1.RowCopy(i,asdfe1);
			blah1.Multx(asdfe1,prode1);
			te1.SetRow(i,prode1);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds2; i++)
		{
			te2.RowCopy(i,asdfe2);
			blah1.Multx(asdfe2,prode2);
			te2.SetRow(i,prode2);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds3; i++)
		{
			te3.RowCopy(i,asdfe3);
			blah1.Multx(asdfe3,prode3);
			te3.SetRow(i,prode3);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds4; i++)
		{
			te4.RowCopy(i,asdfe4);
			blah1.Multx(asdfe4,prode4);
			te4.SetRow(i,prode4);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds5; i++)
		{
			te5.RowCopy(i,asdfe5);
			blah1.Multx(asdfe5,prode5);
			te5.SetRow(i,prode5);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds6; i++)
		{
			te6.RowCopy(i,asdfe6);
			blah1.Multx(asdfe6,prode6);
			te6.SetRow(i,prode6);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds7; i++)
		{
			te7.RowCopy(i,asdfe7);
			blah1.Multx(asdfe7,prode7);
			te7.SetRow(i,prode7);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds8; i++)
		{
			te8.RowCopy(i,asdfe8);
			blah1.Multx(asdfe8,prode8);
			te8.SetRow(i,prode8);
		}
//		fAllEdgeBonds = tedge2;
		fBonds = tedge2;
		fAllBulkBonds = tbulk;
 		fEdgeBonds1 = te1;
		fEdgeBonds2 = te2;
		fEdgeBonds3 = te3;
		fEdgeBonds4 = te4;
		fEdgeBonds5 = te5;
		fEdgeBonds6 = te6;
		fEdgeBonds7 = te7;
		fEdgeBonds8 = te8;		
	}	
	else if (fNormalCode == 5)	// rotate [-1,0,0] to [0,0,-1]
	{
		tedge2 = tempedge;
		tbulk = temp_bulk;
 		te1 = temp_edge1;
		te2 = temp_edge2;
		te3 = temp_edge3;
		te4 = temp_edge4;
		te5 = temp_edge5;
		te6 = temp_edge6;
		te7 = temp_edge7;
		te8 = temp_edge8;	
		blah1 = RotationMatrixB(piby2);
		for (int i = 0; i < kEAMFCC3DEdgeBonds; i++)
		{
			tedge2.RowCopy(i,asdfe);	// take bond
			blah1.Multx(asdfe,prode);		// rotate bond via rotation matrix
			tedge2.SetRow(i,prode);	// place new bond back into temp_bonds
		}
		for (int i = 0; i < kEAMFCC3DBulkBonds; i++)
		{
			tbulk.RowCopy(i,asdfb);
			blah1.Multx(asdfb,prodb);
			tbulk.SetRow(i,prodb);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds1; i++)
		{
			te1.RowCopy(i,asdfe1);
			blah1.Multx(asdfe1,prode1);
			te1.SetRow(i,prode1);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds2; i++)
		{
			te2.RowCopy(i,asdfe2);
			blah1.Multx(asdfe2,prode2);
			te2.SetRow(i,prode2);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds3; i++)
		{
			te3.RowCopy(i,asdfe3);
			blah1.Multx(asdfe3,prode3);
			te3.SetRow(i,prode3);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds4; i++)
		{
			te4.RowCopy(i,asdfe4);
			blah1.Multx(asdfe4,prode4);
			te4.SetRow(i,prode4);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds5; i++)
		{
			te5.RowCopy(i,asdfe5);
			blah1.Multx(asdfe5,prode5);
			te5.SetRow(i,prode5);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds6; i++)
		{
			te6.RowCopy(i,asdfe6);
			blah1.Multx(asdfe6,prode6);
			te6.SetRow(i,prode6);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds7; i++)
		{
			te7.RowCopy(i,asdfe7);
			blah1.Multx(asdfe7,prode7);
			te7.SetRow(i,prode7);
		}
		for (int i = 0; i < kEAMFCC3DEdgeBonds8; i++)
		{
			te8.RowCopy(i,asdfe8);
			blah1.Multx(asdfe8,prode8);
			te8.SetRow(i,prode8);
		}
//		fAllEdgeBonds = tedge2;
		fBonds = tedge2;
		fAllBulkBonds = tbulk;
 		fEdgeBonds1 = te1;
		fEdgeBonds2 = te2;
		fEdgeBonds3 = te3;
		fEdgeBonds4 = te4;
		fEdgeBonds5 = te5;
		fEdgeBonds6 = te6;
		fEdgeBonds7 = te7;
		fEdgeBonds8 = te8;		
	}	

	/* scale to correct lattice parameter */				     		
//	fAllEdgeBonds *= fLatticeParameter;
	fBonds *= fLatticeParameter;
	fAllBulkBonds *= fLatticeParameter;
 	fEdgeBonds1 *= fLatticeParameter;
	fEdgeBonds2 *= fLatticeParameter;
	fEdgeBonds3 *= fLatticeParameter;
	fEdgeBonds4 *= fLatticeParameter;
	fEdgeBonds5 *= fLatticeParameter;
	fEdgeBonds6 *= fLatticeParameter;
	fEdgeBonds7 *= fLatticeParameter;
	fEdgeBonds8 *= fLatticeParameter;	
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
