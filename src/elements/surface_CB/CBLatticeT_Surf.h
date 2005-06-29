/* $Id: CBLatticeT_Surf.h,v 1.1 2005-06-29 22:38:16 hspark Exp $ */
/* created: paklein (12/02/1996)*/
#ifndef _EAMLATTICET_H_
#define _EAMLATTICET_H_

/* base class */
#include "BondLatticeT_Surf.h"

/* direct members */
#include "iArrayT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "dMatrixT.h"

namespace Tahoe {

/** bond vector information needed for Cauchy-Born calculations */
class CBLatticeT_Surf: public BondLatticeT_Surf
{
public:

	/** constructor */
	CBLatticeT_Surf(void);

	/** fetch bond component tensor (R_I R_J R_K R_L) in reduced index
	 * form */
	void BondComponentTensor4(int numbond, dMatrixT& matrix) const;

	/** fetch bond component tensor (R_I R_J) */
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

} // namespace Tahoe 
#endif /* _EAMLATTICET_H_ */
