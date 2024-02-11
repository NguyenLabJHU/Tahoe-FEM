/* $Id: PMLMatT.h,v 1.5 2004-07-15 08:28:02 paklein Exp $ */
/* created: TDN (5/31/2001) */
#ifndef _PML_H_
#define _PML_H_
 
#include "SolidMaterialT.h"
#include "IsotropicT.h"
//#include "Material2DT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

/* forward declarations */
class PMLT;
class C1FunctionT;

class ifstreamT;

/** base class for small strain linear elastic viscoelastic 
 * constitutive law */
class PMLMatT: public SolidMaterialT, IsotropicT//, Material2DT
{
	public:

	/* constructor */
	PMLMatT(ifstreamT& in, const MaterialSupportT& support, const PMLT& element);

	/* apply pre-conditions at the current time step */
	void InitStep(void);

	/* initialize history variable */
	bool NeedsPointInitialization(void) const {return true;}; // declare true
	void PointInitialize(void);                // assigns storage space
	
	/* update/reset internal variables */
	void UpdateHistory(void); // element at a time
	void ResetHistory(void);  // element at a time
	void Load(ElementCardT& element, int ip);
	void Store(ElementCardT& element, int ip);
					
	/* spatial description */
	const dMatrixT& c_ijkl(void); // spatial tangent moduli
	const dSymMatrixT& s_ij(void); // Cauchy stress

	/* material description */
	const dMatrixT& C_IJKL(void); // material tangent moduli
	const dSymMatrixT& S_IJ(void); // PK2 stress
			
	 protected:
	 
	/*calculates relaxation time*/
	virtual const double& DampFacta(dArrayT& ip_coords);	
	virtual const double& DampFactb(dArrayT& ip_coords);
		
	/* returns the strain energy density for the specified strain */
	virtual double StrainEnergyDensity(void);
	
	/*number of spatial dimension*/
	int fNumSD;
	int fNumDOF;
		
	/*stress*/
	dSymMatrixT fStress;
	/*moduli*/
	dMatrixT fModulus;  	
	
	/*Internal state variables*/	
	/*fh - denotes the overstress while 
	 *fs - is the inelastic stress, defined as the the modulus 
	 *of the Maxwell element times the total strain*/

	/*history variables*/
	dSymMatrixT   fStressa_n1;
	dSymMatrixT   fStressa_n;
	dSymMatrixT   fStressb_n1;
	dSymMatrixT   fStressb_n;
	dSymMatrixT   fStress0a_n1;
	dSymMatrixT   fStress0b_n1;
	dSymMatrixT   fStress0a_n;
	dSymMatrixT   fStress0b_n;

	/*number of state variables*/
	int fnstatev;
	
	/* Internal state variables array*/
	dArrayT fstatev;

  	/*Time increment*/
 	const double& fdt; 

	enum FunType {klinear = 1, kquadratic = 2};
 	
 	const PMLT& fPMLElement;
 	
 	C1FunctionT* fFuna;
 	C1FunctionT* fFunb;
 	 		
	double fDampa;
	double fDampb;
	
	double fRefCoorda;
	double fRefCoordb;
	
  };

} // namespace Tahoe 
#endif /*_PML_H_*/
