/* $Id: CBLatticeT.h,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (12/02/1996)                                          */
/* CBLatticeT.h                                                           */

#ifndef _EAMLATTICET_H_
#define _EAMLATTICET_H_

/* base class */
#include "BondLatticeT.h"

/* direct members */
#include "iArrayT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "dMatrixT.h"

class CBLatticeT: public BondLatticeT
{
public:

	/* constructor */
	CBLatticeT(int numlatticedim, int numspatialdim, int numbonds);
	
	/* the Q matrix passed into this constructor is used to rotate the
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
	CBLatticeT(const dMatrixT& Q, int numspatialdim, int numbonds);
			
	/* fetch bond component tensor (R_I R_J R_K R_L) in reduced index
	 * form */
	void BondComponentTensor4(int numbond, dMatrixT& matrix) const;

	/* fetch bond component tensor (R_I R_J) */
	void BondComponentTensor2(int numbond, dArrayT& vector) const;
	void BatchBondComponentTensor2(dArray2DT& comptable) const;
	  		
private:

	/* building the bond component tensors */
	void BondTensor4_2D(const dArrayT& comps, dMatrixT& matrix) const;	
	void BondTensor4_3D(const dArrayT& comps, dMatrixT& matrix) const;	

	void BondTensor2_2D(const dArrayT& comps, dArrayT& vector) const;	
	void BondTensor2_3D(const dArrayT& comps, dArrayT& vector) const;
	
	/* batched versions */	
	void BatchBondTensor2_2D(dArray2DT& comptable) const;	
	void BatchBondTensor2_3D(dArray2DT& comptable) const;
};

#endif /* _EAMLATTICET_H_ */
