/* $Id: SSJ2LinHardT.h,v 1.4 2003-08-08 22:56:06 thao Exp $ */
/* created: paklein (02/12/1997)                                          */
/* Interface for a elastoplastic material that is linearly                */
/* isotropically elastic subject to the Huber-von Mises yield             */
/* condition as fYield with kinematic/isotropic hardening laws            */
/* given by:                                                              */
/* 		H(a) = (1 - ftheta) fH_bar a                                         */
/* K(a) = fYield + ftheta fH_bar a                                        */
/* 		where a is the internal hardening variable                           */
/* 	Note: all calculations are peformed in 3D.                            */

#ifndef _SS_J2_LIN_HARD_T_H_
#define _SS_J2_LIN_HARD_T_H_

/* base class */
/* direct members */
#include "SSJ2LinHardBaseT.h"

namespace Tahoe {

class SSJ2LinHardT: public SSJ2LinHardBaseT
{
public:

	/* constructor */
	SSJ2LinHardT(ifstreamT& in, const SSMatSupportT& support);

	/*free energy density*/
	virtual double StrainEnergyDensity(void);
	virtual double Pressure(void) const;

	/*evaluates stress and consistent modulus*/
	virtual const dMatrixT& c_ijkl(void);
	virtual const dSymMatrixT& s_ij(void);
 
       /*return internal dissipation variables*/
        virtual const dArrayT& InternalStressVars(void);
        virtual const dArrayT& InternalStrainVars(void);
    
	/*material outputs*/
	virtual int NumOutputVariables(void) const;
	virtual void OutputLabels(ArrayT<StringT>& labels) const;
	virtual void ComputeOutput(dArrayT& output);
			
protected:

    /*stres and modulus*/
    dSymMatrixT fStress;
    dMatrixT fModulus;    
};

} // namespace Tahoe 
#endif /* _SS_J2_LIN_HARD_T_H_ */
