/* $Id: WLC.h,v 1.3 2005-02-25 18:56:32 paklein Exp $ */
/* created: TDN (01/22/2001) */

#ifndef _WLC_
#define _WLC_

#include "FSSolidMatT.h"

namespace Tahoe {

class WLC: public FSSolidMatT
{
   public:
  
	/* constructor/destructor */
	WLC(void);

	/* print parameters */
//	virtual void Print(ostream& out) const;
	
	/* strain energy density */
	virtual double StrainEnergyDensity(void);
	virtual const dMatrixT& c_ijkl(void);
	virtual const dSymMatrixT& s_ij(void);
	virtual const dMatrixT& C_IJKL(void);
	virtual const dSymMatrixT& S_IJ(void);
	
	virtual double Pressure(void) const {
		cout << "\n WLC::Pressure: not implemented" << endl;
		throw ExceptionT::kGeneralFail;
		return 0.0;
	};

	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);

	protected:
	/*compute deformed bond lengths*/
	void Compute_R(dArrayT& r, const dSymMatrixT& C);
	void Compute_eigs(dArrayT& eigs, const dSymMatrixT& C);
	
   protected:
	/* return values */
	dMatrixT	fModulus;
	dSymMatrixT     fStress;

   private:  
	/*material parameters*/
	double fN;
	double fA;
	double fT;
	
	/*unit cell parameters*/
	double fl1;
	double fl2;
	double fl3;
	dArrayT f_e1;
	dArrayT f_e2;
	dArrayT f_e3;

	double fR0; /*undeformed RMS distance*/
	double fL;  /*contour length*/
	
	/*bulk parameters*/
	double fgamma;
	double fbeta;
	
	/*work spaces*/
	dSymMatrixT fC;
	dSymMatrixT fStretch;
	dSymMatrixT fStress3D;
	dMatrixT fModulus3D;
	dArrayT fR;
	dArrayT fEigs;
	dArrayT fP;
	ArrayT<dSymMatrixT> fM;
};
}
#endif /* _WLC_ */
