/* $Id: SSJ2LinHard2D.h,v 1.1 2003-05-12 23:38:36 thao Exp $ */
/* created: paklein (02/12/1997)                                          */
/* Interface for a elastoplastic material that is linearly                */
/* isotropically elastic subject to the Huber-von Mises yield             */
/* condition as fYield with kinematic/isotropic hardening laws            */
/* given by:                                                              */
/* 		H(a) = (1 - ftheta) fH_bar a                                         */
/* K(a) = fYield + ftheta fH_bar a                                        */
/* 		where a is the internal hardening variable                           */
/* 	Note: all calculations are peformed in 3D.                            */

#ifndef _SS_J2_LIN_HARD_2D_H_
#define _SS_J2_LIN_HARD_2D_H_

/* base class */
/* direct members */
#include "SSJ2LinHardBaseT.h"
#include "Material2DT.h"

namespace Tahoe {

class SSJ2LinHard2D: public SSJ2LinHardBaseT, public Material2DT
{
public:

	/* constructor */
	SSJ2LinHard2D(ifstreamT& in, const SSMatSupportT& support);

	virtual double Pressure(void) const {
	  cout << "\nSSJ2LinHardT::Pressure: not implemented" <<endl;
	  throw ExceptionT::kGeneralFail;
	  return 0.0;
	};

	/*free energy density*/
	virtual double StrainEnergyDensity(void);

	/*evaluates stress and consistent modulus*/
	virtual const dMatrixT& c_ijkl(void);
	virtual const dSymMatrixT& s_ij(void);
    
    /*material outputs*/
    virtual int NumOutputVariables(void) const;
    virtual void OutputLabels(ArrayT<StringT>& labels) const;
    virtual void ComputeOutput(dArrayT& output);
			
protected:

    /*stress and modulus*/
    dSymMatrixT fStress;
    dMatrixT fModulus;  
    
    /*work space*/
    dSymMatrixT fStrain3D;
    dSymMatrixT fStress3D;
    dMatrixT fModulus3D;  
};

} // namespace Tahoe 
#endif /* _SS_J2_LIN_HARD_2DT_H_ */
