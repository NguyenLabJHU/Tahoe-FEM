/* $Id: SSSimoViscoT.h,v 1.6 2004-07-15 08:28:46 paklein Exp $ */
/* created: TDN (5/31/2001) */
#ifndef _SS_SIMO_VISCO_H_
#define _SS_SIMO_VISCO_H_
 
#include "SSSolidMatT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

/*small strain linear viscoelastic constitutive law */
class SSSimoViscoT: public SSSolidMatT
{
	public:

	/*constructor*/
	SSSimoViscoT(ifstreamT& in, const SSMatSupportT& support);

	/** return the pressure associated with the last call to 
	 * SolidMaterialT::s_ij. \note NOT IMPLEMENTED */
	virtual double Pressure(void) const {
		ExceptionT::GeneralFail("SSSimoViscoT::Pressure", "not implemented");
		return 0.0;
	};

	/* apply pre-conditions at the current time step */
	virtual void InitStep(void){SSSolidMatT::InitStep();}

	/*initialize history variable*/
	virtual bool NeedsPointInitialization(void) const {return true;}; // declare true
	virtual void PointInitialize(void);                // assigns storage space
	
	/* update/reset internal variables */
	virtual void UpdateHistory(void); // element at a time
	virtual void ResetHistory(void);  // element at a time
	void Load(ElementCardT& element, int ip);
	void Store(ElementCardT& element, int ip);

	/*inquire if dissipation variables used in material force calculation are needed*/
	virtual bool HasDissipVar(void) const {return true;}
	virtual const iArrayT& InternalDOF(void) const;
			
	protected:
	
	enum Spring {kEquilibrium = 0, kNonEquilibrium};
	protected:
	 
	/*Internal state variables*/	
	/*fh - denotes the overstress while 
	 *fs - is the inelastic stress, defined as the the modulus 
	 *of the Maxwell element times the total strain*/

	/*preceding values*/		 
	dSymMatrixT   fdevQ_n;
	dSymMatrixT   fdevSin_n;
	dArrayT        fmeanQ_n;
	dArrayT        fmeanSin_n;
	
	/*current values*/
	dSymMatrixT   fdevQ;
	dSymMatrixT   fdevSin;
	dArrayT   fmeanQ;
	dArrayT   fmeanSin;

	/*internal dissipation variables*/
	iArrayT fInternalDOF;
	dSymMatrixT fViscStrain;
	dSymMatrixT fViscStress;

	/* Internal state variables array*/
	int fnstatev;
	dArrayT fstatev;

 	/*relaxation times*/
	double ftauS;
	double ftauB;

	/* dt/tau*/
	double fndtS;
	double fndtB;	
  };

 inline const iArrayT& SSSimoViscoT::InternalDOF(void) const 
 {
	return fInternalDOF;
 }
 
} // namespace Tahoe 
#endif /*_SS_SIMO_VISCO_H_*/
