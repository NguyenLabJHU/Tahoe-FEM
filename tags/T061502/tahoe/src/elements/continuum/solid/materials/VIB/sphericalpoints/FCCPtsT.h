/* $Id: FCCPtsT.h,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (03/26/1999)                                          */
/* FCC lattice of points                                                  */

#ifndef _FCC_PTS_T_H_
#define _FCC_PTS_T_H_

/* base class */
#include "SpherePointsT.h"

class FCCPtsT: public SpherePointsT
{
public:

	/* constructor */
	FCCPtsT(int num_shells, double bond_length);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;	

	/* generate sphere points:
	 *
	 *   phi   = angle about z from x
	 *   theta = angle about x from z
	 *
	 * The final orientation is generated by applied the
	 * phi and theta rotations in succession about the local
	 * axes.
	 */
	virtual const dArray2DT& SpherePoints(double phi, double theta);

private:

	/* set lattice coordinates */
	void SetCoords(void);

private:

	/* parameters */
	int	   fNumShells;
	double fBondLength;	// nearest neighbor bond length		
};

#endif /* _FCC_PTS_T_H_ */