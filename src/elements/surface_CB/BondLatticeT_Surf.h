/* $Id: BondLatticeT_Surf.h,v 1.1 2005-06-29 22:38:16 hspark Exp $ */
/* created: paklein (01/07/1997) */
#ifndef _BondLatticeT_Surf_H_
#define _BondLatticeT_Surf_H_

/* direct members */
#include "iArrayT.h"
#include "dArrayT.h"
#include "dSymMatrixT.h"
#include "dArray2DT.h"
#include "dMatrixT.h"

namespace Tahoe {

/** container for bond information */
class BondLatticeT_Surf
{
public:

	/** constructor must be followed by call to BondLatticeT_Surf::Initialize to
	 * initialize the bond table information */
	BondLatticeT_Surf(void);
	
	/** destructor */
	virtual ~BondLatticeT_Surf(void);

	/** The Q matrix is used to rotate the
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
	void Initialize(const dMatrixT* Q = NULL);

	/** \name accessors */
	/*@{*/
	const iArrayT& BondCounts(void) const;
	const dArrayT& DeformedLengths(void) const;
	dArrayT& DeformedLengths(void);
	const dArray2DT& Bonds(void) const;
//	int NumberOfLatticeDim(void) const;
//	int NumberOfSpatialDim(void) const;
	int NumberOfBonds(void) const { return fBonds.MajorDim(); };
	dSymMatrixT& Stretch(void) { return fStretch; };
	/*@}*/

	/* compute deformed bond lengths from the given Green strain */
	void ComputeDeformedLengths(const dSymMatrixT& strain);

protected:

	/** initialize bond table values */
	virtual void LoadBondTable(void) = 0;
	
protected:

	iArrayT		fBondCounts;
	dArray2DT	fBonds;			/* undeformed bond vector components */
	dArrayT 	fDefLength;		/* list of deformed bond lengths */
	dMatrixT	fQ;				/* bond vector transformation matrix */

	/** \name work space */
	/*@{*/
	dArrayT		fBondSh;		/**< shallow bond vector */
	dArrayT 	fBondDp;		/**< deep bond vector */
//	dMatrixT	fLatDimMatrix;	/**< matrix with same dimensions as lattice */
	dSymMatrixT	fStrain;		/**< needed if LatticeDim != SpatialDim */  		
	dSymMatrixT	fStretch;		/**< stretch tensor */
	/*@}*/
};

/* inlines */
inline const iArrayT& BondLatticeT_Surf::BondCounts(void) const { return fBondCounts; }
inline const dArrayT& BondLatticeT_Surf::DeformedLengths(void) const { return fDefLength; }
inline dArrayT& BondLatticeT_Surf::DeformedLengths(void) { return fDefLength; }
inline const dArray2DT& BondLatticeT_Surf::Bonds(void) const { return fBonds; }

} /* namespace Tahoe */

#endif /* _BondLatticeT_Surf_H_ */
