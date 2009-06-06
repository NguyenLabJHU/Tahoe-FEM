/* $Id: BondLatticeT.h,v 1.10 2009-06-06 17:28:49 hspark Exp $ */
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
	const iArrayT& EdgeBulkCounts(void) const;
	const iArrayT& Edge1Counts(void) const;
	const iArrayT& Edge2Counts(void) const;
	const iArrayT& Edge3Counts(void) const;
	const iArrayT& Edge4Counts(void) const;
	const iArrayT& Edge5Counts(void) const;
	const iArrayT& Edge6Counts(void) const;
	const iArrayT& Edge7Counts(void) const;
	const iArrayT& Edge8Counts(void) const;	
	
	/* Surf3 arrays for {110} surfaces */
	const iArrayT& Surf3Counts(void) const;
	const dArrayT& DeformedLengths(void) const;
	const dArrayT& DeformedSurfaceLengths(void) const;
	const dArrayT& DeformedEdgeLengths(void) const;
	const dArrayT& DeformedBulk(void) const;
	const dArrayT& DeformedSurf1(void) const;
	const dArrayT& DeformedSurf2(void) const;
	const dArrayT& DeformedSurf3(void) const;
	const dArrayT& DeformedEdge1(void) const;
	const dArrayT& DeformedEdge2(void) const;
	const dArrayT& DeformedEdge3(void) const;
	const dArrayT& DeformedEdge4(void) const;
	const dArrayT& DeformedEdge5(void) const;
	const dArrayT& DeformedEdge6(void) const;
	const dArrayT& DeformedEdge7(void) const;
	const dArrayT& DeformedEdge8(void) const;	
	const dArrayT& DeformedEdgeBulk(void) const;
	dArrayT& DeformedLengths(void);
	dArrayT& DeformedSurfaceLengths(void);
	dArrayT& DeformedEdgeLengths(void);
	const dArray2DT& Bonds(void) const;
//	int NumberOfLatticeDim(void) const;
//	int NumberOfSpatialDim(void) const;
	int NumberOfBonds(void) const { return fBonds.MajorDim(); };
	dSymMatrixT& Stretch(void) { return fStretch; };
	iArrayT& AtomTypes(void) {return fAtomType; };
	iArrayT& EdgeTypes(void) {return fEdgeType; };	
	/*@}*/

	/* compute deformed bond lengths from the given Green strain */
	void ComputeDeformedLengths(const dSymMatrixT& strain);
	void ComputeDeformedSurfaceLengths(const dSymMatrixT& strain);
	void ComputeDeformedEdgeLengths(const dSymMatrixT& strain);

	/* Similar deformed lengths function but for representative bulk/surface atoms */
	void ComputeDeformedBulkBonds(const dSymMatrixT& strain);
	void ComputeDeformedSurf1Bonds(const dSymMatrixT& strain);
	void ComputeDeformedSurf2Bonds(const dSymMatrixT& strain);
	void ComputeDeformedSurf3Bonds(const dSymMatrixT& strain);

	/* EDGE STUFF */
	void ComputeDeformedEdgeBulkBonds(const dSymMatrixT& strain);
	void ComputeDeformedEdge1Bonds(const dSymMatrixT& strain);
	void ComputeDeformedEdge2Bonds(const dSymMatrixT& strain);
	void ComputeDeformedEdge3Bonds(const dSymMatrixT& strain);
	void ComputeDeformedEdge4Bonds(const dSymMatrixT& strain);
	void ComputeDeformedEdge5Bonds(const dSymMatrixT& strain);
	void ComputeDeformedEdge6Bonds(const dSymMatrixT& strain);
	void ComputeDeformedEdge7Bonds(const dSymMatrixT& strain);
	void ComputeDeformedEdge8Bonds(const dSymMatrixT& strain);

protected:

	/** initialize bond table values */
	virtual void LoadBondTable(void) = 0;
	
protected:

	iArrayT		fBondCounts;
	iArrayT     fBulkCounts;
	iArrayT     fSurf1Counts;
	iArrayT     fSurf2Counts;
	iArrayT 	fSurf3Counts;
	iArrayT		fBulkBondCounts;	// EDGE
	iArrayT		fEdge1Counts;
	iArrayT		fEdge2Counts;
	iArrayT		fEdge3Counts;
	iArrayT		fEdge4Counts;
	iArrayT		fEdge5Counts;
	iArrayT		fEdge6Counts;
	iArrayT		fEdge7Counts;
	iArrayT		fEdge8Counts;	
	dArray2DT	fBonds;			/* undeformed bond vector components */
	dArray2DT   fBulkBonds;		/* undeformed bond lengths for a representative bulk atom */
	dArray2DT   fSurf1Bonds;    /* undeformed bond lengths for a representative surface atom */
	dArray2DT   fSurf2Bonds;    /* undeformed bond lengths for a representative 1 layer into the bulk atom */
	dArray2DT 	fSurf3Bonds;	/* undeformed bond lengths for a representative 2 layer into the bulk atom */
	dArray2DT  	fEdgeBonds1;	/* undeformed bond lengths for type-1 edge atom */
	dArray2DT  	fEdgeBonds2;	/* undeformed bond lengths for type-2 edge atom */
	dArray2DT  	fEdgeBonds3;	/* undeformed bond lengths for type-3 edge atom */
	dArray2DT  	fEdgeBonds4;	/* undeformed bond lengths for type-4 edge atom */
	dArray2DT  	fEdgeBonds5;	/* undeformed bond lengths for type-5 edge atom */
	dArray2DT  	fEdgeBonds6;	/* undeformed bond lengths for type-6 edge atom */
	dArray2DT  	fEdgeBonds7;	/* undeformed bond lengths for type-7 edge atom */
	dArray2DT  	fEdgeBonds8;	/* undeformed bond lengths for type-8 edge atom */		
	dArray2DT	fAllBulkBonds;	/* undeformed bond lengths for a representative bulk atom - edge case */
	dArrayT 	fDefLength;		/* list of deformed bond lengths */
	dArrayT     fDefBulk;		/* list of deformed bulk bonds */
	dArrayT     fDefSurf1;		/* list of deformed surface bonds */
	dArrayT     fDefSurf2;      /* list of deformed bonds for atom 1 layer into the bulk */
	dArrayT		fDefSurf3;		/* list of deformed bonds for atom 2 layers into the bulk */
	dArrayT		fDefBulkLength;	/* list of deformed bulk bonds (edge case) */
	dArrayT		fDefEdge1;		/* list of deformed type-1 edge atom bond lengths */
	dArrayT		fDefEdge2;		/* list of deformed type-2 edge atom bond lengths */
	dArrayT		fDefEdge3;		/* list of deformed type-3 edge atom bond lengths */
	dArrayT		fDefEdge4;		/* list of deformed type-4 edge atom bond lengths */
	dArrayT		fDefEdge5;		/* list of deformed type-5 edge atom bond lengths */
	dArrayT		fDefEdge6;		/* list of deformed type-6 edge atom bond lengths */
	dArrayT		fDefEdge7;		/* list of deformed type-7 edge atom bond lengths */
	dArrayT		fDefEdge8;		/* list of deformed type-8 edge atom bond lengths */		
	dMatrixT	fQ;				/* bond vector transformation matrix */
	iArrayT     fAtomType;		/* interaction indicator type (0-6) for surface CB */
	/* 0=s1/s1, 1=s1/s2, 2=s1/bulk, 3=s2/s1, 4=s2/s2, 5=s2/bulk, 6=bulk/bulk */
	iArrayT 	fEdgeType;		/* interaction indicator type for edge CB */
			
	/** \name work space */
	/*@{*/
	dArrayT		fBondSh;		/**< shallow bond vector */
	dArrayT     fBondShB;		/**< shallow bond vector for bulk atom */
	dArrayT     fBondShS;		/**< shallow bond vector for all surface atoms */
	dArrayT     fBondShS1;		/**< shallow bond vector for surface atom 1 */
	dArrayT     fBondShS2;		/**< shallow bond vector for surface atom 2 */
	dArrayT		fBondShS3;		/**< shallow bond vector for surface atom 3 */
	dArrayT		fBondShBE;		/**< shallow bond vector for bulk atom for edge CB case */
	dArrayT		fBondShE;		/**< shallow bond vector for all edge atoms */
	dArrayT		fBondShE1;		/**< shallow bond vector for type-1 edge atoms */	
	dArrayT		fBondShE2;		/**< shallow bond vector for type-2 edge atoms */
	dArrayT		fBondShE3;		/**< shallow bond vector for type-3 edge atoms */
	dArrayT		fBondShE4;		/**< shallow bond vector for type-4 edge atoms */
	dArrayT		fBondShE5;		/**< shallow bond vector for type-5 edge atoms */
	dArrayT		fBondShE6;		/**< shallow bond vector for type-6 edge atoms */
	dArrayT		fBondShE7;		/**< shallow bond vector for type-7 edge atoms */
	dArrayT		fBondShE8;		/**< shallow bond vector for type-8 edge atoms */	
	dArrayT 	fBondDp;		/**< deep bond vector */
	dArrayT     fBondDpS;		/**< deep bond vector for all surface atoms */
	dArrayT     fBondDpB;		/**< deep bond vector for bulk atom */
	dArrayT     fBondDpS1;		/**< deep bond vector for surface atom 1 */
	dArrayT     fBondDpS2;		/**< deep bond vector for surface atom 2 */
	dArrayT		fBondDpS3;		/**< deep bond vector for surface atom 3 */
	dArrayT		fBondDpBE;		/**< deep bond vector for bulk atom (edge CB case) */
	dArrayT		fBondDpE;		/**< deep bond vector for all edge atoms */
	dArrayT		fBondDpE1;		/**< deep bond vector for type-1 edge atoms */
	dArrayT		fBondDpE2;		/**< deep bond vector for type-2 edge atoms */
	dArrayT		fBondDpE3;		/**< deep bond vector for type-3 edge atoms */
	dArrayT		fBondDpE4;		/**< deep bond vector for type-4 edge atoms */
	dArrayT		fBondDpE5;		/**< deep bond vector for type-5 edge atoms */
	dArrayT		fBondDpE6;		/**< deep bond vector for type-6 edge atoms */
	dArrayT		fBondDpE7;		/**< deep bond vector for type-7 edge atoms */
	dArrayT		fBondDpE8;		/**< deep bond vector for type-8 edge atoms */		
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
inline const iArrayT& BondLatticeT::EdgeBulkCounts(void) const { return fBulkBondCounts; }
inline const iArrayT& BondLatticeT::Edge1Counts(void) const { return fEdge1Counts; }
inline const iArrayT& BondLatticeT::Edge2Counts(void) const { return fEdge2Counts; }
inline const iArrayT& BondLatticeT::Edge3Counts(void) const { return fEdge3Counts; }
inline const iArrayT& BondLatticeT::Edge4Counts(void) const { return fEdge4Counts; }
inline const iArrayT& BondLatticeT::Edge5Counts(void) const { return fEdge5Counts; }
inline const iArrayT& BondLatticeT::Edge6Counts(void) const { return fEdge6Counts; }
inline const iArrayT& BondLatticeT::Edge7Counts(void) const { return fEdge7Counts; }
inline const iArrayT& BondLatticeT::Edge8Counts(void) const { return fEdge8Counts; }
inline const dArrayT& BondLatticeT::DeformedLengths(void) const { return fDefLength; }
inline dArrayT& BondLatticeT::DeformedLengths(void) { return fDefLength; }
inline const dArrayT& BondLatticeT::DeformedSurfaceLengths(void) const { return fDefLength; }
inline dArrayT& BondLatticeT::DeformedSurfaceLengths(void) {return fDefLength; }
inline const dArrayT& BondLatticeT::DeformedEdgeLengths(void) const { return fDefLength; }
inline dArrayT& BondLatticeT::DeformedEdgeLengths(void) {return fDefLength; }
inline const dArrayT& BondLatticeT::DeformedBulk(void) const { return fDefBulk; }
inline const dArrayT& BondLatticeT::DeformedSurf1(void) const { return fDefSurf1; }
inline const dArrayT& BondLatticeT::DeformedSurf2(void) const { return fDefSurf2; }
inline const dArrayT& BondLatticeT::DeformedSurf3(void) const { return fDefSurf3; }
inline const dArrayT& BondLatticeT::DeformedEdgeBulk(void) const { return fDefBulkLength; }
inline const dArrayT& BondLatticeT::DeformedEdge1(void) const { return fDefEdge1; }
inline const dArrayT& BondLatticeT::DeformedEdge2(void) const { return fDefEdge2; }
inline const dArrayT& BondLatticeT::DeformedEdge3(void) const { return fDefEdge3; }
inline const dArrayT& BondLatticeT::DeformedEdge4(void) const { return fDefEdge4; }
inline const dArrayT& BondLatticeT::DeformedEdge5(void) const { return fDefEdge5; }
inline const dArrayT& BondLatticeT::DeformedEdge6(void) const { return fDefEdge6; }
inline const dArrayT& BondLatticeT::DeformedEdge7(void) const { return fDefEdge7; }
inline const dArrayT& BondLatticeT::DeformedEdge8(void) const { return fDefEdge8; }
inline const dArray2DT& BondLatticeT::Bonds(void) const { return fBonds; }

} /* namespace Tahoe */

#endif /* _BONDLATTICET_H_ */
