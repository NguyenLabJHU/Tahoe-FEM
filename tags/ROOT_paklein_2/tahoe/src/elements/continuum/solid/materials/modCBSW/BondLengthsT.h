/* $Id: BondLengthsT.h,v 1.2 2002-07-02 19:56:06 cjkimme Exp $ */
/* created: paklein (05/20/1997)                                          */
/* Class to compute/manage all bond lengths and derivatives               */
/* for the 2 unit cell, diamond cubic, modified Cauchy-Born,              */
/* constitutive equations.                                                */

#ifndef _BOND_LENGTHS_T_H_
#define _BOND_LENGTHS_T_H_

/* direct members */
#include "dArrayT.h"
#include "dMatrixT.h"
#include "dArray2DT.h"
#include "dSymMatrixT.h"


namespace Tahoe {

class BondLengthsT
{
public:

	/* Constructor */
	BondLengthsT(const dMatrixT& Q); //pass empty Q for default orientation

	/* Destructor */
	virtual ~BondLengthsT(void);

	/* set free dof - triggers recomputation */
	void SetdXsi(const dMatrixT& CIJ, const dArrayT& Xsi); //computes R + Xsi
	void SetdC(const dMatrixT& CIJ);
	void SetAll(const dMatrixT& CIJ);

	/* Deformed lengths */
	const dArrayT& Lengths(void) const;

	/* Derivatives wrt. Xsi */
	const dArray2DT& dl_dXsi(void) const;
	const dArray2DT& d2l_dXsidXsi(void) const;

	/* Derivatives wrt. C */
	const dArray2DT& dl_hat_dC(void) const;
	const dArray2DT& d2l_hat_dCdC(void) const;
	const dArray2DT& d2l_hat_dCdXsi(void) const;

protected:

	/* Compute dR/dC and d2R/dCdC */
	void dR(const dArray2DT& R, const dArrayT& l, const dMatrixT& C,
		dArray2DT& dC, dArray2DT& dCdC);

private:

	/* form symmetrized contribution to d2l_dCdXsi */
	void Symmetrized(dMatrixT& mat, const dArrayT& Rmod);

	/* called by constructor */
	void Initialize(const dMatrixT& Q);

protected:
	  	    	  	
	/* Undeformed bond vectors */
	dArray2DT	fR;
	dArray2DT	fRmod;
	
	/* Deformed bond lengths */
	dArrayT		fl;
	
	/* Derivatives wrt Xsi */
	dArray2DT	fdl_dXsi;
	dArray2DT	fd2l_dXsidXsi;

	/* Derivatives wrt C */
	dArray2DT	fdl_dC;
	dArray2DT	fd2l_dCdC;

	/* Derivatives mixed wrt Xsi and C */
	dArray2DT	fd2l_dCdXsi;	
	
	/* work space */
	dArrayT		fTempVec;
	dMatrixT	fTempMat1;
	dMatrixT	fTempRank4;
	dSymMatrixT	fTempSymMat1;
	dMatrixT	fTempMixed;
	
};

} // namespace Tahoe 
#endif /* _BOND_LENGTHS_T_H_ */
