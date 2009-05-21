/* $Id: RGSplitT.h,v 1.6 2009-05-21 22:30:27 tdnguye Exp $ */
/* created: TDN (01/22/2001) */
#ifndef _RGSplitT_
#define _RGSplitT_

/**S Reese, S Govindjee (1998) IJSS 35:3455-3482:
 **This is an implementation that will be removed in the near future.  
 **Use the new more general implementation RGSplitT2
 **\

/* base class */
#include "RGViscoelasticityT.h"

namespace Tahoe {

class RGSplitT: public RGViscoelasticityT
{
   public:
  
	/* constructor/destructor */
	RGSplitT(void);

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

	/*free energy density*/
	virtual double Energy(const dArrayT& lambda_bar, const double J, const int type);

	/*calculates principal values of deviatoric Kirchoff stress given principal values of deviatoric stretch tensor*/
	virtual void DevStress(const dArrayT& lambda_bar, dArrayT& tau, const int type);
	/*calculates mean Kirchhoff stress tensor given J*/
	virtual double MeanStress(const double J, const int type);

	/*calculates principal values of deviatoric stiffness given principal values of deviatoric stretch tensor*/
 	virtual void DevMod(const dArrayT& lambda_bar,dSymMatrixT& eigenmodulus, const int type);
	/*calculates bulk mod given J*/
	virtual double MeanMod(const double J, const int type);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
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
	
//   private:  
	/* spectral operations */
	SpectralDecompT fSpectralDecompSpat;

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
	
  	/*viscosities*/
	double fietaS;
	double fietaB;
	
	/*moduli for NeoHookean Potential*/
	double fmu_eq;
	double fmu_neq;
	double fkappa_eq;
	double fkappa_neq;
};
}
#endif /* _RGSplitT_ */
