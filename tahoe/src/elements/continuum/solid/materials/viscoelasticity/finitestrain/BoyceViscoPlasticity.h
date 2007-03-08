/* $Id: BoyceViscoPlasticity.h,v 1.1 2007-03-08 18:13:01 tdnguye Exp $ */
/* created: TDN (01/22/2001) */
#ifndef _BoyceViscoPlasticity_
#define _BoyceViscoPlasticity_

/**Reese and Govindjee IJSS 1998:  nonlinear viscoelasticity model with  **
 **constant isotropic viscosity tensor. Volumetric/Deviatoric split formulation/ **
 **Constitutive relation is calculated by default using a compressible Neo-Hookean 
 **potential.  Other potentials can be implemented in derived classes by overloading 
 **BoyceViscoPlasticity::dWdE and BoyceViscoPlasticity::dWdE.    **\

/* base class */
#include "BoyceBaseT.h"
#include "InvLangevin.h"
namespace Tahoe {

class BoyceViscoPlasticity: public BoyceBaseT
{
   public:
  
	/* constructor/destructor */
	BoyceViscoPlasticity(void);

	/* initialize, update/reset internal variables */
	virtual void PointInitialize(void);              
	virtual void UpdateHistory(void); // element at a time
	virtual void ResetHistory(void);  // element at a time

	/* apply pre-conditions at the current time step */
	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	
	/* strain energy density */
	virtual double StrainEnergyDensity(void);
	virtual const dMatrixT& c_ijkl(void);
	virtual const dSymMatrixT& s_ij(void);
	virtual const dMatrixT& C_IJKL(void);
	virtual const dSymMatrixT& S_IJ(void);

	/*compute output variables*/ 
	virtual int NumOutputVariables() const; 
	virtual void OutputLabels(ArrayT<StringT>& labels) const; 
	virtual void ComputeOutput(dArrayT& output);

   protected:
   virtual void Initialize(void);
   
   private:
	virtual void Compute_Calg(const dArrayT& eigenstretch, const dArrayT& eigenstretch_e);
	virtual void ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e);    
	dMatrixT DefGrad(void);

   protected:
	/* return values */
	dMatrixT fModulus;
	dSymMatrixT fStress;

	/*inverse langevin function*/
	InvLangevin fInvL;
	
   private:  
	/* spectral operations */
	SpectralDecompT fSpectralDecompSpat;

	/*internal variable*/
	double* fs;
	double* fs_n;
	
	/*work space*/
	dMatrixT fF3D;
	dSymMatrixT fStretch;
	dMatrixT fInverse;
	dSymMatrixT fStress3D;

	dMatrixT fFe;
	dSymMatrixT fStretch_e;
	dMatrixT fRe;
	
	dArrayT     fEigs;
	dArrayT     fEigs_e;
	
	dArrayT     fEigs_Stress;

	dArrayT fRes;
	dArrayT fdel;
	
	dMatrixT fK;
	dMatrixT fH;
	dMatrixT fG;
	dMatrixT fM;	
	dMatrixT fCalg;
	dMatrixT    fModulus3D;
	dMatrixT    fModMat;
	
	/*back stress shear modulus*/
 	double fmu_R;
	/*limiting chain stretch*/
	double flambda_L;

 	/*ref strain rate*/
	double fgammadot0;
	/*yield stress*/
	double fs0;
	/*saturation strength*/
	double fs_ss;
	/*hardening stiffness*/
	double fh;
	/*stress stiffening modulus*/
	double fA;
	/*pressure param*/
	double falpha;
	/*temperature*/
	double fT;

	
	/*moduli for NeoHookean Potential*/
	double fmu;
	double fkappa;
};
}
#endif /* _BoyceViscoPlasticity_ */
