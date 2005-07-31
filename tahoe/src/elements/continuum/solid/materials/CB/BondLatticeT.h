/* $Id: BondLatticeT.h,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (01/07/1997)                                          */
/* BondLatticeT.h                                                         */

#ifndef _BONDLATTICET_H_
#define _BONDLATTICET_H_

/* direct members */
#include "iArrayT.h"
#include "dArrayT.h"
#include "dSymMatrixT.h"
#include "dArray2DT.h"
#include "dMatrixT.h"

class BondLatticeT
{
public:

	/* constructor */
	BondLatticeT(int numlatticedim, int numspatialdim, int numbonds);
	
	/* The Q matrix passed into this constructor is used to rotate the
	 * bond vectors into the orientation prescribed by Q.  No check is
	 * performed on the orthogonality of Q, only its dimensions.  Q is
	 * deep copied.  Q is defined as:
	 *
	 *			Q = d x_natural / d x_global
	 *
	 * So that the vectors are transformed by:
	 *
	 *			r_global = Transpose[Q].r_natural
	 */
	BondLatticeT(const dMatrixT& Q, int numspatialdim, int numbonds);
	
	/* destructor */
	virtual ~BondLatticeT(void);

	/* initialize bond table */
	void Initialize(void);
		
	/* references to the lattice lists */
	const iArrayT& BondCounts(void) const;
	const dArrayT& DeformedLengths(void) const;

	/* accessors */
	int NumberOfLatticeDim(void) const;
	int NumberOfSpatialDim(void) const;
	int NumberOfBonds(void) const;

protected:

	/* initialize bond table values */
	virtual void LoadBondTable(void) = 0;
	
	/* compute deformed bond lengths from the given Green strain */
	void ComputeDeformedLengths(const dSymMatrixT& strain);
	
protected:

	int fIsInitialized;
	int fNumLatticeDim;	/* dim of the bond vectors */
	int	fNumSpatialDim;	/* dim of the model geometry */
	int fNumBonds;
	
	iArrayT		fBondCounts;
	dArray2DT	fBonds;			/* undeformed bond vector components */
	dArrayT 	fDefLength;		/* list of deformed bond lengths */
	dMatrixT	fQ;				/* bond vector transformation matrix */

	/* work space */
	dArrayT		fBondSh;		/* shallow bond vector */
	dArrayT 	fBondDp;		/* deep bond vector */
	dMatrixT	fLatDimMatrix;	/* matrix with same dimensions as lattice */
	dSymMatrixT	fStrain;		/* needed if LatticeDim != SpatialDim */  		
};

#endif /* _BONDLATTICET_H_ */
