/* created: Majid T. Manzari (04/16/2003) */
#ifndef _MR_SS_KSTV_H_
#define _MR_SS_KSTV_H_

/* base classes */
#include "SSIsotropicMatT.h"
#include "HookeanMatT.h"
#include "MRSSNLHardT.h"

namespace Tahoe {

class MRSSKStV: public SSIsotropicMatT,
				public HookeanMatT,
				public MRSSNLHardT
{
  public:

	/* constructor */
	MRSSKStV(ifstreamT& in, const SSMatSupportT& support);
	
	/* initialization */
	virtual void Initialize(void);

	/* form of tangent matrix (symmetric by default) */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/* update internal variables */
	virtual void UpdateHistory(void);

	/* reset internal variables to last converged solution */
	virtual void ResetHistory(void);

	/** \name spatial description */
	/*@{*/
	/** spatial tangent modulus */
	virtual const dMatrixT& c_ijkl(void);

	virtual const dMatrixT& cdisc_ijkl(void);

	/** Cauchy stress */
	virtual const dSymMatrixT& s_ij(void);

	/** return the pressure associated with the last call to 
	 * SolidMaterialT::s_ij. See SolidMaterialT::Pressure
	 * for more information. */
	virtual double Pressure(void) const { return fStress.Trace()/3.0; };
	/*@}*/

	/* returns the strain energy density for the specified strain */
	virtual double StrainEnergyDensity(void);

	/* returns the number of variables computed for nodal extrapolation
	 * during for element output, ie. internal variables */
	virtual int  NumOutputVariables(void) const;
	virtual void OutputLabels(ArrayT<StringT>& labels) const;
	virtual void ComputeOutput(dArrayT& output);

        /*
         * Test for localization using "current" values for Cauchy
         * stress and the spatial tangent moduli. Returns 1 if the
         * determinant of the acoustic tensor is negative and returns
         * the normal for which the determinant is minimum. Returns 0
         * of the determinant is positive.
         */
	 int IsLocalized(dArrayT& normal);

protected:

	/* set modulus */
 	virtual void SetModulus(dMatrixT& modulus); 
         int loccheck;
 
  private:
  
  	/* return values */
  	dSymMatrixT	fStress;
  	dMatrixT	fModulus, fModulusCe;
    dMatrixT    fModulusdisc;

};

} // namespace Tahoe 
#endif /* _MR_SS_KSTV_H_ */
