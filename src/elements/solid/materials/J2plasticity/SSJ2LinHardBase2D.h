/* $Id: SSJ2LinHardBase2D.h,v 1.2 2004-07-15 08:28:34 paklein Exp $ */
/* created: paklein (02/12/1997)                                          */
/* Interface for a elastoplastic material that is linearly                */
/* isotropically elastic subject to the Huber-von Mises yield             */
/* condition as fYield with kinematic/isotropic hardening laws            */
/* given by:                                                              */
/* 		H(a) = (1 - ftheta) fH_bar a                                         */
/* K(a) = fYield + ftheta fH_bar a                                        */
/* 		where a is the internal hardening variable                           */
/* 	Note: all calculations are peformed in 3D.                            */

#ifndef _SS_J2_LIN_HARD_BASE_2D_H_
#define _SS_J2_LIN_HARD_BASE_2D_H_

/* base class */
/* direct members */
#include "SSIsotropicMatT.h"
#include "dSymMatrixT.h"
#include "dMatrixT.h"
#include "dArrayT.h"

namespace Tahoe {

/* forward declarations */
class ElementCardT;

class SSJ2LinHardBase2D: public SSIsotropicMatT
{
public:

        /* constructor */
        SSJ2LinHardBase2D(ifstreamT& in, const SSMatSupportT& support);
  
	/* output name */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;
	
	/* apply pre-conditions at the current time step */
	virtual void InitStep(void){SSSolidMatT::InitStep();}
       
	/*manage history variables*/
 	/*initialize history variable*/
	virtual bool NeedsPointInitialization(void) const {return true;}; // declare true
	virtual void PointInitialize(void);
	virtual void UpdateHistory(void);
	virtual void ResetHistory(void);
	virtual void Load(ElementCardT& element, int ip);
	virtual void Store(ElementCardT& element, int ip);

        /*inquire if dissipation variables used in material force calculation are n
eeded*/
        virtual bool HasDissipVar(void) const {return true;}
        virtual const iArrayT& InternalDOF(void) const;

protected:
	/*radial return mapping*/
	double YieldCondition(const dSymMatrixT& relstress, double alpha) const;
	bool PlasticLoading(const dSymMatrixT& devtrialstress);
	/* return the correction to stress vector computed by the mapping the
	 * stress back to the yield surface, if needed */
	const dSymMatrixT& StressCorrection(const dSymMatrixT& devtrialstress);
	
	/* return the correction to moduli due to plasticity (if any)
	 *
	 * Note: Return mapping occurs during the call to StressCorrection.
	 *       The element passed in is already assumed to carry current
	 *       internal variable values */
	const dMatrixT& ModuliCorrection(void);

	/*accessors*/
	double  H(double a) const;
	const double dH(void) const;
	double  K(double a) const;
	const double dK(void) const;
			
protected:

	/*material properties*/
	double fMu;     /*elastic shear stiffness*/
	double fKappa;  /*elastic bulk stiffness*/
	double fYield;	/* initial flow stress (fYield > 0) */
	double fH_bar;	/* hardening parameter (fH_bar > 0) */	
	double ftheta;	/* (0 < ftheta < 1) 				*/

	const double fthird;

	/*flag*/
	bool fplastic;

	/* element level internal variables */
	double fnstatev;             // number of internal state variables 
	dArrayT fstatev;             // internal state variable array

	/*current values*/
	dArrayT falpha;              // equivalent plastic strain
	dSymMatrixT fBeta;           // back stress 
	dSymMatrixT fPlasticStrain;  // deviatoric part of the plastic strain
    
	/*previous values*/
	dArrayT falpha_n;              // equivalent plastic strain
	dSymMatrixT fBeta_n;           // back stress 
	dSymMatrixT fPlasticStrain_n;  // deviatoric part of the plastic strain

       /*internal dissipation variables*/
        dArrayT fInternalStressVars;
        dArrayT fInternalStrainVars;
        iArrayT fInternalDOF;
        
	/* work space */
	dSymMatrixT fUnitNorm;      //unit normal to the stress surface
	double ftrial;
	double fdgamma;
	dSymMatrixT	fStressCorr;
	dMatrixT	fModuliCorr;
		
	dSymMatrixT fRelStress;
	dSymMatrixT fDevStrain; /* deviatoric part of the strain tensor */
	dMatrixT   fTensorTemp; 
};

 inline const iArrayT& SSJ2LinHardBase2D::InternalDOF(void) const
 {
   return fInternalDOF;
 }

 inline double SSJ2LinHardBase2D::H(double a) const
 {
   return ( (1.0 - ftheta)*fH_bar*a );
 }

 inline const double SSJ2LinHardBase2D::dH(void) const
 {
   return ( (1.0 - ftheta)*fH_bar );
 }

 inline double SSJ2LinHardBase2D::K(double a) const
 {
   return ( fYield + ftheta*fH_bar*a );
 }

 inline const double SSJ2LinHardBase2D::dK(void) const
 {
   return ( ftheta*fH_bar );
 }

} // namespace Tahoe 
#endif /* _SS_J2_LIN_HARD_Base_2D_H_ */
