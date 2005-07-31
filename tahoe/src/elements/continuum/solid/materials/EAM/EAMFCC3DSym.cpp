/* $Id: EAMFCC3DSym.cpp,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (12/06/1996)                                          */
/* EAMFCC3DSym.cpp                                                        */

#include "EAMFCC3DSym.h"

/* constructor */
EAMFCC3DSym::EAMFCC3DSym(ifstreamT& in, int EAMcode, int numspatialdim, int numbonds):
	EAMFCC3D(in, EAMcode, numspatialdim, numbonds)
{

}

EAMFCC3DSym::EAMFCC3DSym(ifstreamT& in, const dMatrixT& Q, int EAMcode, int numspatialdim,
	int numbonds):
	EAMFCC3D(in, Q, EAMcode, numspatialdim, numbonds)
{

}

/**********************************************************************
* Protected
**********************************************************************/
	
void EAMFCC3DSym::LoadBondTable(void)
{
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

	/* dimension check */				     		
	if (fBonds.MajorDim() != fNumBonds ||
	    fBonds.MinorDim() != kEAMFCC3DNumLatticeDim) throw eGeneralFail;
	
	/* copy into reference table */
	for (int i = 0; i < fNumBonds; i++)
		for (int j = 0; j < kEAMFCC3DNumLatticeDim; j++)
			fBonds(i,j) = bonddata[i][j];
			
	/* scale to correct lattice parameter */				     		
	fBonds *= fLatticeParameter;
}
