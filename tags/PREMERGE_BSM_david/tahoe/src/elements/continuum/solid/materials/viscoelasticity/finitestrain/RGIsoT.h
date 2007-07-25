/* $Id: RGIsoT.h,v 1.2 2006-08-21 16:46:24 tdnguye Exp $ */
/* created: TDN (01/22/2001) */
#ifndef _RGIsoT_
#define _RGIsoT_

/**Reese and Govindjee IJSS 1998:  nonlinear viscoelasticity model with  **
 **constant isotropic viscosity tensor;  no volumetric/deviatoric split  **
 **Derived types need only to overload RGIsoT::dWdE and RGIsoT::dWdE.    **\
 
 
/* base class */
#include "RGViscoelasticityT.h"

/* direct members */
#include "SpectralDecompT.h"

namespace Tahoe {

class RGIsoT: public RGViscoelasticityT
{
   public:
  
	/* constructor/destructor */
	RGIsoT(void);
	
	/*compute stress and moduli*/
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
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

   protected:

   	/* principal values of the PK2 stress given principal values of the total and elastic stretch 
	 * tensors, i.e., the principal stretches squared 
	 * process index refers the nonequilibrium process number.  process_index = -1 denotes eq. process.*/
	virtual void dWdE(const dArrayT& eigenstretch, dArrayT& eigenstress, const int process_index) = 0;

	/*calculates the derivative 1/lambda_B pdf(S_A)(lambda_B)*/
	virtual void ddWddE(const dArrayT& eigenstretch, dArrayT& eigenstress,
		 dSymMatrixT& eigenmod, const int process_index) = 0;

	/* return true of model is purely 2D, plain stress */
	virtual bool PurePlaneStress(void) const { return false; };

   private:
	void ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e, const int process_index);
	void ComputeiKAB(dSymMatrixT& eigenmodulus_neq);
    
   protected:

	/* return values */
	dMatrixT fModulus;
	dSymMatrixT fStress;
	
   private:  
	/* spectral operations */
	SpectralDecompT fSpectralDecompSpat;
	
	/*work space*/
	dSymMatrixT fbe;
	dSymMatrixT fb_tr;
	dSymMatrixT fb;
	dMatrixT	fF3D;
	
	dArrayT     fEigs;
	dArrayT     fEigs_e;
	dArrayT     fEigs_tr;
	
	dArrayT		fdWdE_eq;
	dArrayT		fdWdE_neq;
	dSymMatrixT	fddWddE_eq;	
	dSymMatrixT	fddWddE_neq;

	dSymMatrixT    fStress3D;
	dSymMatrixT    fInverse;
	dMatrixT	fCalg;
	dMatrixT    fModulus3D;
	dMatrixT    fModMat;
  	dMatrixT    fiKAB;
	
  	/*viscosities*/
	double fietaS;
	double fietaB;
};
}
#endif /* _RGIsoT_ */
