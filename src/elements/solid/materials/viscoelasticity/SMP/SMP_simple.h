/* $Id: SMP_simple.h,v 1.3 2006-11-14 22:58:06 thao Exp $ */
/* created: TDN (01/22/2001) */
#ifndef _SMP_simple_
#define _SMP_simple_

/* base class */
#include "RGSplitT.h"

namespace Tahoe {

class SMP_simple: public RGSplitT
{
   public:

	enum EnergyType {kMooneyRivlin=0,
					kLangevin=1};
					 
	enum ViscType {kNone = -1, 
		       kSimple=0, 
					kPower=1,
					kBergStromBoyce=2}; 
  
	/* constructor/destructor */
	SMP_simple(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;
	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	
	virtual void InitStep(void);

	/*compute output variables*/ 
	virtual int NumOutputVariables() const; 
	virtual void OutputLabels(ArrayT<StringT>& labels) const; 
	virtual void ComputeOutput(dArrayT& output);

	/*viscosity*/
	virtual void Viscosity(double& ietaS, double& ietaB, const int type);
	
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

   private:
	/* set inverse of thermal transformation - return true if active */
//	virtual bool SetInverseThermalTransformation(dMatrixT& F_trans_inv);  			

	virtual void Compute_Calg(const dArrayT& tau_dev, const dSymMatrixT& dtau_dev, const double& tau_m, 
						const double& dtau_m, dMatrixT& Calg, const int type);

	virtual void ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e, 
	                   dArrayT& eigenstress, dSymMatrixT& eigenmodulus, const int process_num);
    
   protected:
	/*Reference Temperature*/
	double fRefTemperature;
	double fTg;
//	dArrayT fTemperature;
	
	/*moduli for Mooney Rivlin Potential, n_process x 2 (c1, c2)*/
    /*rubber moduli at reference temperature
	double fc1_eq;         
	double fc2_eq;         
	double fgamma_eq;
	*/
	dArray2DT fPot;
	dArray2DT fVisc;
	int fPotType;
	int fViscType;
};
}
#endif /* _SMP_simple_ */
