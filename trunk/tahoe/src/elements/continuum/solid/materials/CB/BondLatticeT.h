/* $Id: BondLatticeT.h,v 1.9 2007-06-12 16:03:17 hspark Exp $ */
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

/** container for bond information */
class BondLatticeT
{
public:

	/** constructor must be followed by call to BondLatticeT::Initialize to
	 * initialize the bond table information */
	BondLatticeT(void);
	
	/** destructor */
	virtual ~BondLatticeT(void);

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
	const iArrayT& BulkCounts(void) const;
	const iArrayT& Surf1Counts(void) const;
	const iArrayT& Surf2Counts(void) const;
	
	/* Surf3 arrays for {110} surfaces */
	const iArrayT& Surf3Counts(void) const;
	const dArrayT& DeformedLengths(void) const;
	const dArrayT& DeformedBulk(void) const;
	const dArrayT& DeformedSurf1(void) const;
	const dArrayT& DeformedSurf2(void) const;
	const dArrayT& DeformedSurf3(void) const;
	dArrayT& DeformedLengths(void);
	const dArray2DT& Bonds(void) const;
//	int NumberOfLatticeDim(void) const;
//	int NumberOfSpatialDim(void) const;
	int NumberOfBonds(void) const { return fBonds.MajorDim(); };
	dSymMatrixT& Stretch(void) { return fStretch; };
	iArrayT& AtomTypes(void) {return fAtomType; };
	/*@}*/

	/* compute deformed bond lengths from the given Green strain */
	void ComputeDeformedLengths(const dSymMatrixT& strain);

	/* Similar deformed lengths function but for representative bulk/surface atoms */
	void ComputeDeformedBulkBonds(const dSymMatrixT& strain);
	void ComputeDeformedSurf1Bonds(const dSymMatrixT& strain);
	void ComputeDeformedSurf2Bonds(const dSymMatrixT& strain);
	void ComputeDeformedSurf3Bonds(const dSymMatrixT& strain);

protected:

	/** initialize bond table values */
	virtual void LoadBondTable(void) = 0;
	
protected:

	iArrayT		fBondCounts;
	iArrayT     fBulkCounts;
	iArrayT     fSurf1Counts;
	iArrayT     fSurf2Counts;
	iArrayT 	fSurf3Counts;
	dArray2DT	fBonds;			/* undeformed bond vector components */
	dArray2DT   fBulkBonds;		/* undeformed bond lengths for a representative bulk atom */
	dArray2DT   fSurf1Bonds;    /* undeformed bond lengths for a representative surface atom */
	dArray2DT   fSurf2Bonds;    /* undeformed bond lengths for a representative 1 layer into the bulk atom */
	dArray2DT 	fSurf3Bonds;	/* undeformed bond lengths for a representative 2 layer into the bulk atom */
	dArrayT 	fDefLength;		/* list of deformed bond lengths */
	dArrayT     fDefBulk;		/* list of deformed bulk bonds */
	dArrayT     fDefSurf1;		/* list of deformed surface bonds */
	dArrayT     fDefSurf2;      /* list of deformed bonds for atom 1 layer into the bulk */
	dArrayT		fDefSurf3;		/* list of deformed bonds for atom 2 layers into the bulk */
	dMatrixT	fQ;				/* bond vector transformation matrix */
	iArrayT     fAtomType;		/* interaction indicator type (0-6) for surface CB */
	/* 0=s1/s1, 1=s1/s2, 2=s1/bulk, 3=s2/s1, 4=s2/s2, 5=s2/bulk, 6=bulk/bulk */
			
	/** \name work space */
	/*@{*/
	dArrayT		fBondSh;		/**< shallow bond vector */
	dArrayT     fBondShB;		/**< shallow bond vector for bulk atom */
	dArrayT     fBondShS1;		/**< shallow bond vector for surface atom 1 */
	dArrayT     fBondShS2;		/**< shallow bond vector for surface atom 2 */
	dArrayT		fBondShS3;		/**< shallow bond vector for surface atom 3 */
	dArrayT 	fBondDp;		/**< deep bond vector */
	dArrayT     fBondDpB;		/**< deep bond vector for bulk atom */
	dArrayT     fBondDpS1;		/**< deep bond vector for surface atom 1 */
	dArrayT     fBondDpS2;		/**< deep bond vector for surface atom 2 */
	dArrayT		fBondDpS3;		/**< deep bond vector for surface atom 3 */
//	dMatrixT	fLatDimMatrix;	/**< matrix with same dimensions as lattice */
	dSymMatrixT	fStrain;		/**< needed if LatticeDim != SpatialDim */  		
	dSymMatrixT	fStretch;		/**< stretch tensor */
	/*@}*/
};

/* inlines */
inline const iArrayT& BondLatticeT::BondCounts(void) const { return fBondCounts; }
inline const iArrayT& BondLatticeT::BulkCounts(void) const { return fBulkCounts; }
inline const iArrayT& BondLatticeT::Surf1Counts(void) const { return fSurf1Counts; }
inline const iArrayT& BondLatticeT::Surf2Counts(void) const { return fSurf2Counts; }
inline const iArrayT& BondLatticeT::Surf3Counts(void) const { return fSurf3Counts; }
inline const dArrayT& BondLatticeT::DeformedLengths(void) const { return fDefLength; }
inline dArrayT& BondLatticeT::DeformedLengths(void) { return fDefLength; }
inline const dArrayT& BondLatticeT::DeformedBulk(void) const { return fDefBulk; }
inline const dArrayT& BondLatticeT::DeformedSurf1(void) const { return fDefSurf1; }
inline const dArrayT& BondLatticeT::DeformedSurf2(void) const { return fDefSurf2; }
inline const dArrayT& BondLatticeT::DeformedSurf3(void) const { return fDefSurf3; }
inline const dArray2DT& BondLatticeT::Bonds(void) const { return fBonds; }

} /* namespace Tahoe */

#endif /* _BONDLATTICET_H_ */
