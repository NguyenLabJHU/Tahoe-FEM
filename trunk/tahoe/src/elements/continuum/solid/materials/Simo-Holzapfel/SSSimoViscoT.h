/* $Id: SSSimoViscoT.h,v 1.1 2002-10-05 02:48:19 thao Exp $ */
/* created:   TDN (5/31/2001) */
#ifndef _SS_SIMO_VISCO_H_
#define _SS_SIMO_VISCO_H_
 
#include "SSStructMatT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

class ifstreamT;

/*small strain linear viscoelastic constitutive law */
class SSSimoViscoT: public SSStructMatT
{
	public:

	/*constructor*/
	SSSimoViscoT(ifstreamT& in, const SmallStrainT& element);

	/** return the pressure associated with the last call to 
	 * StructuralMaterialT::s_ij. \note NOT IMPLEMENTED */
	virtual double Pressure(void) const {
		cout << "\n SSSimoViscoT::Pressure: not implemented" << endl;
		throw eGeneralFail;
		return 0.0;
	};
	
	/*print parameters*/
	void Print(ostream& out) const;
	void PrintName(ostream& out) const;
		
	/* apply pre-conditions at the current time step */
	void InitStep(void){SSStructMatT::InitStep();}

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
			
	protected:
	
	enum Spring {kEquilibrium = 0, kNonEquilibrium};
	protected:

	/*Form Residual flag*/
	const GlobalT::StateT& fRunState;
	 
	/*Internal state variables*/	
	/*fh - denotes the overstress while 
	 *fs - is the inelastic stress, defined as the the modulus 
	 *of the Maxwell element times the total strain*/

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

	/*number of state variables*/
	int fnstatev;

	/* Internal state variables array*/
	dArrayT fstatev;

  	/*Time increment*/
 	const double& fdt; 
	
	/*modulus*/
	dMatrixT fModulus;
	/*stress*/
	dSymMatrixT fStress;

 	 /*relaxation times*/
	double ftauS;
	double ftauB;

	/* dt/tau*/
	double fndtS;
	double fndtB;	
  };

} // namespace Tahoe 
#endif /*_SS_SIMO_VISCO_H_*/
