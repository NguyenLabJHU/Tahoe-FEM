/* $Id: PMLMatT.h,v 1.1 2002-01-21 06:49:02 thao Exp $ */
/* created:   TDN (5/31/2001) */

#ifndef _PML_H_
#define _PML_H_
 
#include "StructuralMaterialT.h"
#include "PMLT.h"
#include "IsotropicT.h"
#include "Material2DT.h"
#include "C1FunctionT.h"
class ifstreamT;

/** base class for small strain linear elastic viscoelastic 
 * constitutive law */
class PMLMatT: public StructuralMaterialT,IsotropicT,Material2DT
{
	public:

	/*constructor*/
	PMLMatT(ifstreamT& in, const PMLT& element);
	
	/*print parameters*/
	void Print(ostream& out) const;
	void PrintName(ostream& out) const;
		
	/*calculate instantaneous moduli*/
	virtual void Initialize(void);
		
	/* apply pre-conditions at the current time step */
	void InitStep(void);

	/*initialize history variable*/
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
	
	/*Form Residual flag*/
	const GlobalT::StateT& fRunState;

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

#endif /*_PML_H_*/
