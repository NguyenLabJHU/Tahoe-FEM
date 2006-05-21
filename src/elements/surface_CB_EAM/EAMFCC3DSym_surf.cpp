/* $Id: EAMFCC3DSym_surf.cpp,v 1.1 2006-05-21 15:55:19 hspark Exp $ */
/* created: paklein (12/06/1996) */
#include "EAMFCC3DSym_surf.h"

/* bond parameters */
const int kEAMFCC3DSymNumBonds = 27;
const int kEAMFCC3DNumLatticeDim =  3;

using namespace Tahoe;

/* constructor */
EAMFCC3DSym_surf::EAMFCC3DSym_surf(void) {
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

	/* all bonds appear twice */
	fBondCounts = 2;
	
	/* clear deformed lengths for now */
	fDefLength = 0.0;

	/*
	 * Undeformed bond data for unit cube - the table includes
	 * all neighbors in a (2 x 2 x 2) unit cell area except those
	 * located at [±1 ±1 ±1], ie the body diagonal neighbors.
	 */
	double bonddata[kEAMFCC3DSymNumBonds][kEAMFCC3DNumLatticeDim] =
	
						   {/* face centers - 1st nearest */
				/*  1 */	{ 0.5, 0.5, 0.0},
				/*  2 */	{-0.5, 0.5, 0.0},
				/*  3 */	{ 0.5, 0.0, 0.5},
				/*  4 */	{-0.5, 0.0, 0.5},
				/*  5 */	{ 0.0, 0.5, 0.5},
				/*  6 */	{ 0.0,-0.5, 0.5},
	  			
	  			/* face corners - 4th nearest */
				/*  7 */	{ 1.0, 1.0, 0.0},
				/*  8 */	{-1.0, 1.0, 0.0},
				/*  9 */	{ 1.0, 0.0, 1.0},
				/* 10 */	{-1.0, 0.0, 1.0},
				/* 11 */	{ 0.0, 1.0, 1.0},
				/* 12 */	{ 0.0,-1.0, 1.0},
				
							/* edge corners - 2nd nearest */
				/* 13 */	{ 1.0, 0.0, 0.0},
				/* 14 */	{ 0.0, 1.0, 0.0},
				/* 15 */	{ 0.0, 0.0, 1.0},
					
							/* far face corners - 3rd nearest */					
				/* 16 */	{ 1.0, 0.5, 0.5},
				/* 17 */	{ 1.0, 0.5,-0.5},
				/* 18 */	{ 1.0,-0.5, 0.5},
				/* 19 */	{ 1.0,-0.5,-0.5},
				/* 20 */	{ 0.5, 1.0, 0.5},
				/* 21 */	{ 0.5, 1.0,-0.5},
				/* 22 */	{-0.5, 1.0, 0.5},
				/* 23 */	{-0.5, 1.0,-0.5},
				/* 24 */	{ 0.5, 0.5, 1.0},
				/* 25 */	{ 0.5,-0.5, 1.0},
				/* 26 */	{-0.5, 0.5, 1.0},
				/* 27 */	{-0.5,-0.5, 1.0}
				     		};
	
	/* copy into reference table */
	for (int i = 0; i < kEAMFCC3DSymNumBonds; i++)
		for (int j = 0; j < 3; j++)
			fBonds(i,j) = bonddata[i][j];
			
	/* scale to correct lattice parameter */				     		
	fBonds *= fLatticeParameter;
}
