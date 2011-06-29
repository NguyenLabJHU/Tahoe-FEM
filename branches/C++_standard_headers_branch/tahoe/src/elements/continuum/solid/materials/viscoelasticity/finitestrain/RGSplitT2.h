/* $Id: RGSplitT2.h,v 1.3 2009-05-21 22:30:27 tdnguye Exp $ */
/* created: TDN (01/22/2001) */
#ifndef _RGSplitT2_
#define _RGSplitT2_

/**S Reese, S Govindjee (1998) IJSS 35:3455-3482:  
 **nonlinear viscoelasticity model with  **
 **constant isotropic viscosity tensor. Volumetric/Deviatoric split formulation/ **
 **Constitutive relation is calculated by default using a compressible Neo-Hookean 
 **potential.  Other potentials can be implemented in derived classes by overloading 
 **RGSplitT2::dWdE and RGSplitT2::dWdE.    **\

/* base class */
#include "RGViscoelasticityT.h"
#include "PotentialT.h"
#include "C1FunctionT.h"

namespace Tahoe {

class RGSplitT2: public RGViscoelasticityT
{
   public:
  
	/* constructor/destructor */
	RGSplitT2(void);

//	enum EnergyType {kEQ=0, kNEQ=1}; 
	
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

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;
	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);

	/**compute mechanical strains*/
	virtual const dMatrixT& MechanicalDeformation(void);

	/**compute thermal strains*/
	virtual const dMatrixT& ThermalDeformation_Inverse(void);
	
	/*@}*/

   protected:
   virtual void Initialize(void);
   
   private:
	virtual void Compute_Calg(const dArrayT& tau_dev, const dSymMatrixT& dtau_dev, const double& tau_m, 
						const double& dtau_m, dMatrixT& Calg, const int type);
	virtual void ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e, 
	                   dArrayT& eigenstress, dSymMatrixT& eigenmodulus, const int process_num);    
   protected:

	/* return values */
	dMatrixT fModulus;
	dSymMatrixT fStress;
	
	/* spectral operations */
	SpectralDecompT fSpectralDecompSpat;

	/*mechanical strains*/
	dMatrixT fF_M;
	
	/*thermal strains*/
	dMatrixT fF_T_inv;
	
	/*work space*/
	dSymMatrixT fb;
	dSymMatrixT fbe;
	dSymMatrixT fb_tr;
	dMatrixT fF3D;
	dSymMatrixT fInverse;
	
	dArrayT     fEigs;
	dArrayT     fEigs_e;
	dArrayT     fEigs_tr;
	dArrayT     fEigs_dev;

	dArrayT	    ftau_EQ;
	dArrayT     ftau_NEQ;

	dSymMatrixT fStress3D;
	dSymMatrixT fDtauDe_EQ;
	dSymMatrixT fDtauDe_NEQ;
	dMatrixT fCalg;

	dMatrixT    fModulus3D;
	dMatrixT    fModMat;
  	dMatrixT    fiKAB;
	
	/*potential*/
	ArrayT<PotentialT*> fPot;
  	/*viscosities*/
	ArrayT<C1FunctionT*> fVisc_s;
	ArrayT<C1FunctionT*> fVisc_b;	
	double fietaS;
	double fietaB;
};

}
#endif /* _RGSplitT2_ */
