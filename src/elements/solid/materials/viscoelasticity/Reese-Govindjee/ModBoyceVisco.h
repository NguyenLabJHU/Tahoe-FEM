/* $Id: ModBoyceVisco.h,v 1.1 2007-07-25 06:43:26 thao Exp $ */
/* created: TDN (01/22/2001) */
#ifndef _ModBoyceVisco_
#define _ModBoyceVisco_

/* base class */
#include "RGSplitT2.h"
#include "InvLangevin.h"

namespace Tahoe {

class ModBoyceVisco: public RGSplitT2
{
   public:
  
	/* constructor/destructor */
	ModBoyceVisco(void);

	/* initialize, update/reset internal variables */
	virtual void PointInitialize(void);              
	virtual void UpdateHistory(void); // element at a time
	virtual void ResetHistory(void);  // element at a time


	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;
	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	
	virtual void InitStep(void);

	/*compute output variables*/ 
	virtual int NumOutputVariables() const; 
	virtual void OutputLabels(ArrayT<StringT>& labels) const; 
	virtual void ComputeOutput(dArrayT& output);

	virtual double ShearViscosity(const double smag, const double sy);	
	

	protected:
	virtual void Initialize(void);

   private:
	/* set inverse of thermal transformation - return true if active */
//	virtual bool SetInverseThermalTransformation(dMatrixT& F_trans_inv);  			

	virtual void Compute_Calg(const dArrayT& tau_dev, const dSymMatrixT& dtau_dev, const double& tau_m, 
						const double& dtau_m, dMatrixT& Calg, const int type);

	virtual void ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e, 
	                   dArrayT& eigenstress, dSymMatrixT& eigenmodulus, const int process_num);
    
   protected:
	
	
	
	/*low temp viscous flow*/
	double fTemperature; /*temperature*/
	double fetaS0;		/*reference shear viscosity*/
	double fQS;			/*activation energy for the viscoplastic flow*/
	double fsy0;		/*initial yield strength*/;
	double fsinf;		/*steady state limit of the yield strength*/
	double fh;			/*hardening modulus*/
		
	/*accessors for internal variables*/
	double* fsy;
	double* fsy_n;
	
	/*Residual*/
	dArrayT fRes;
	dArrayT fDelta;
	dMatrixT fGAB;
	dMatrixT fDAB;
	dMatrixT fDABbar;
	dMatrixT fMat;
};
}
#endif /* _ModBoyceVisco_ */
