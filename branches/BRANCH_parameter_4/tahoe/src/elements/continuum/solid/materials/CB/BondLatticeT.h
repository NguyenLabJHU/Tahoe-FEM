/* $Id: BondLatticeT.h,v 1.4 2004-06-26 05:58:47 paklein Exp $ */
/* created: paklein (01/07/1997) */
#ifndef _BONDLATTICET_H_
#define _BONDLATTICET_H_

/* direct members */
#include "iArrayT.h"
#include "dArrayT.h"
#include "dSymMatrixT.h"
#include "dArray2DT.h"
#include "dMatrixT.h"

namespace Tahoe {

/** container for unit cell bonds */
class BondLatticeT
{
public:

	/** constructor */
	BondLatticeT(int numlatticedim, int numspatialdim, int numbonds);
	
	/** The Q matrix passed into this constructor is used to rotate the
	 * bond vectors into the orientation prescribed by Q.  No check is
	 * performed on the orthogonality of Q, only its dimensions.  Q is
	 * deep copied.  Q is defined as:
	 \f[
	 	\mathbf{Q} = \frac{\partial \mathbf{x}_{natural}}{\partial \mathbf{x}_{global}}
	 \f]
	 * So that the vectors are transformed by:
	 \f[
	 	\mathbf{r}_{global} = \mathbf{Q}^T \mathbf{r}_{natural}
	 \f]
	 */
	BondLatticeT(const dMatrixT& Q, int numspatialdim, int numbonds);
	
	/** destructor */
	virtual ~BondLatticeT(void) {};

	/** initialize bond table */
	void Initialize(void);
		
	/** \name accessors */
	/*@{*/
	const iArrayT& BondCounts(void) const;
	const dArrayT& DeformedLengths(void) const;
	const dArray2DT& Bonds(void) const;
	int NumberOfLatticeDim(void) const;
	int NumberOfSpatialDim(void) const;
	int NumberOfBonds(void) const;
	/*@}*/

	/** compute deformed bond lengths from the given Green strain */
	void ComputeDeformedLengths(const dSymMatrixT& strain);

protected:

	/** initialize bond table values */
	virtual void LoadBondTable(void) = 0;
	
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

inline const iArrayT& BondLatticeT::BondCounts(void) const { return fBondCounts; }
inline const dArrayT& BondLatticeT::DeformedLengths(void) const { return fDefLength; }
inline const dArray2DT& BondLatticeT::Bonds(void) const { return fBonds; }
inline int BondLatticeT::NumberOfLatticeDim(void) const { return fNumLatticeDim; }
inline int BondLatticeT::NumberOfSpatialDim(void) const { return fNumSpatialDim; }
inline int BondLatticeT::NumberOfBonds(void) const { return fNumBonds; }

} /* namespace Tahoe */

#endif /* _BONDLATTICET_H_ */
