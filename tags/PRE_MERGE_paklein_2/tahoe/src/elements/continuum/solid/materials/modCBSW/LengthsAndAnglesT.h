/* $Id: LengthsAndAnglesT.h,v 1.3 2002-07-05 22:28:22 paklein Exp $ */
/* created: paklein (05/26/1997)                                          */
/* Class to compute/manage all bond angles and derivatives                */
/* for the 2 unit cell, diamond cubic, modified Cauchy-Born,              */
/* constitutive equations.                                                */

#ifndef _LENGTHSANDANGLES_T_H_
#define _LENGTHSANDANGLES_T_H_

/* base class */
#include "BondLengthsT.h"

namespace Tahoe {

/* forward declarations */
class iArray2DT;

class LengthsAndAnglesT: public BondLengthsT
{
public:

	/* constructor */
	LengthsAndAnglesT(const dMatrixT& Q, const iArray2DT& pairs);

	/* set free dof - triggers recomputation */
	void SetdXsi(const dMatrixT& CIJ, const dArrayT& Xsi);
	void SetdC(const dMatrixT& CIJ);
	void SetAll(const dMatrixT& CIJ);

	/* Accessors */
	const dArrayT& Cosines(void) const;

	/* Derivatives wrt. Xsi */
	const dArray2DT& dCos_dXsi(void) const;
	const dArray2DT& d2Cos_dXsidXsi(void) const;

	/* Derivatives wrt. C */
	const dArray2DT& dCos_hat_dC(void) const;
	const dArray2DT& d2Cos_hat_dCdC(void) const;
	const dArray2DT& d2Cos_hat_dCdXsi(void) const;

private:

	/* called by constructor */
	void Initialize(void);
	  	
protected:

	/* bondpairs */
	const iArray2DT& fPairs;

	/* Undeformed bond vectors */
	dArray2DT	fR_3;
	
	/* Deformed lengths and cosine */
	dArrayT		fl_3;

	dArray2DT	fdl_3_dC;			// Derivatives wrt C
	dArray2DT	fd2l_3_dCdC;

//Cos(theta) = alpha/beta
	
	/* deformed cosine */
	dArrayT		fCos12;

	/* Derivatives wrt Xsi */
	dArray2DT	fdCos_dXsi;
	dArray2DT	fd2Cos_dXsidXsi;

	/* Derivatives wrt C */
	dArray2DT	fdCos_dC;
	dArray2DT	fd2Cos_dCdC;

	/* Derivatives mixed wrt Xsi and C */
	dArray2DT	fd2Cos_dCdXsi;	

private:

//alpha
	dArrayT		fda_dXsi;		//Derivatives wrt. Xsi
	dMatrixT	fd2a_dXsidXsi;
	dMatrixT	fdah_dC;		//Derivatives wrt. C			
	dMatrixT	fd2ah_dCdC;
	dMatrixT	fd2ah_dCdXsi;	//Derivatives mixed wrt. Xsi and C

//beta
	dArrayT		fdb_dXsi;		//Derivatives wrt. Xsi
	dMatrixT	fd2b_dXsidXsi;		
	dMatrixT	fdbh_dC;		//Derivatives wrt. C
	dMatrixT	fd2bh_dCdC;
	dMatrixT	fd2bh_dCdXsi;	//Derivatives mixed wrt. Xsi and C

	/* work space */
	dMatrixT	fTempMat2;
	dSymMatrixT	fTempSymMat2;
};

} // namespace Tahoe 
#endif /* _LENGTHSANDANGLES_T_H_ */
