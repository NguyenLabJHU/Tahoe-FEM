/* $Id: FDSimoViscoBaseT.h,v 1.1 2003-03-19 19:03:18 thao Exp $ */
/* created:   TDN (5/31/2001) */

#ifndef _FD_SIMO_VISCO_BASE_H_
#define _FD_SIMO_VISCO_BASE_H_
 
#include "FSSolidMatT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

/*small strain linear viscoelastic constitutive law */
class FDSimoViscoBaseT: public FSSolidMatT
{
	public:

	/*constructor*/
	FDSimoViscoBaseT(ifstreamT& in, const FSMatSupportT& support);

	/** return the pressure associated with the last call to 
	 * SolidMaterialT::s_ij. \note NOT IMPLEMENTED */
	virtual double Pressure(void) const {
		cout << "\n FDSimoViscoBaseT::Pressure: not implemented" << endl;
		throw ExceptionT::kGeneralFail;
		return 0.0;
	};

	/*print parameters*/
	void Print(ostream& out) const;
	void PrintName(ostream& out) const;
		
	/* apply pre-conditions at the current time step */
	void InitStep(void){FSSolidMatT::InitStep();};

	/*initialize history variable*/
	bool NeedsPointInitialization(void) const {return true;}; // declare true
	void PointInitialize(void);                // assigns storage space
	
	/* update/reset internal variables */
	void UpdateHistory(void); // element at a time
	void ResetHistory(void);  // element at a time
	void Load(ElementCardT& element, int ip);
	void Store(ElementCardT& element, int ip);

	//compute output variables
	virtual int NumOutputVariables() const;
	virtual void OutputLabels(ArrayT<StringT>& labels) const;
	virtual void ComputeOutput(dArrayT& output);

	virtual double StrainEnergyDensity(void);

	protected:

	enum Spring {kEquilibrium = 0, kNonEquilibrium = 1};
        protected:

	/*number of state variables*/
	int fnstatev;
	
	/* Internal state variables array*/
	dArrayT fstatev;

	/*preceding values*/		 
	/*deviatoric*/
	dSymMatrixT   fDevOverStress_n;
	dSymMatrixT   fDevInStress_n;

	/*mean*/
	dArrayT        fVolOverStress_n;
	dArrayT        fVolInStress_n;
	
	/*current values*/
	/*deviatoric*/
	dSymMatrixT   fDevOverStress;
	dSymMatrixT   fDevInStress;
	/*mean*/
	dArrayT   fVolOverStress;
	dArrayT   fVolInStress;

};
}
#endif /*_FD_SIMO_VISCO_BASE_H_*/
