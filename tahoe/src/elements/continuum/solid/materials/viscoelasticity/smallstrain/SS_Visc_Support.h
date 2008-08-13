/* $Id: SS_Visc_Support.h,v 1.1 2008-08-13 12:45:41 tdnguye Exp $ */
/* created: TDN (5/31/2001) */
#ifndef _SS_VISC_SUPPORT_H_
#define _SS_VISC_SUPPORT_H_
 
#include "dSymMatrixT.h"
#include "ElementCardT.h"
namespace Tahoe {

/** small strain linear viscoelastic constitutive law */
class SS_Visc_Support
{
	public:

	/** constructor */
	SS_Visc_Support(void);

	/* update/reset internal variables */
	virtual void SetSateVariables();
	virtual void Update(ElementCardT& element,  int nip); // element at a time
	virtual void Reset(ElementCardT& element,  int nip);  // element at a time
	void Load(ElementCardT& element, int ip);
	void Store(ElementCardT& element, int ip);
	

	protected:
	void ComputeStress(const dSymMatrixT& strain,  dSymMatrixT& stress, double dt, int index);
	void ComputeStress2D(const dSymMatrixT& strain,  dSymMatrixT& stress, double dt, 
							int constraint, int index);
	
	void SetModulus(dMatrixT& modulus, double dt, int index);
	void SetModulus2D(dMatrixT& modulus, double dt, int constraint, int index);
	
	enum Spring {kEquilibrium = 0, kNonEquilibrium};
	protected:
	 
	/*Internal state variables*/	
	/*fh - denotes the overstress while 
	 *fs - is the inelastic stress, defined as the the modulus 
	 *of the Maxwell element times the total strain*/

	/*preceding values*/		 
	ArrayT<dSymMatrixT>    fdevQ_n;
	ArrayT<dSymMatrixT>    fdevSin_n;
	ArrayT<dArrayT>        fmeanQ_n;
	ArrayT<dArrayT>        fmeanSin_n;
	
	/*current values*/
	ArrayT<dSymMatrixT>   fdevQ;
	ArrayT<dSymMatrixT>   fdevSin;
	ArrayT<dArrayT>       fmeanQ;
	ArrayT<dArrayT>       fmeanSin;

	/* Internal state variables array*/
	int fnumprocess;
	int fnstatev;
	dArrayT fstatev;

	/*strain energy potentials*/ 
	dArrayT fMu;
	dArrayT fKappa;

 	/*relaxation times*/
	dArrayT ftauS;
	dArrayT ftauB;

	private:
	/*workspace*/
	dSymMatrixT fStress3D;
	dSymMatrixT fStrain3D;
  };

} // namespace Tahoe 
#endif /*_SS_VISCO_H_*/
