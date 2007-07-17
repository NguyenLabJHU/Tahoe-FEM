/* $Id: SMP_simple.h,v 1.4 2007-07-17 20:14:00 thao Exp $ */
/* created: TDN (01/22/2001) */
#ifndef _SMP_simple_
#define _SMP_simple_

/* base class */
#include "RGSplitT2.h"
#include "InvLangevin.h"

namespace Tahoe {

class SMP_simple: public RGSplitT2
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

	/*fictive temperature*/
	virtual double FictiveTemperature(const double deltaneq);
	
	/*viscosity*/
	virtual double RetardationTime(const double Temperature, const double deltaneq);
	virtual double ShearViscosity(const double Temperature, const double deltaneq, const double smag, const double sy);	
	
	/**compute thermal strains*/
	virtual const dMatrixT& ThermalDeformation_Inverse(void);

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
	/*Reference Temperature*/
	double fT0; /*Temperature at To*/
	
	/*thermal expansion*/
	double fT2;			/*zero entropy temperature*/
	double falphar;		/*high temp CTE*/
	double falphag;		/*low temp CTE*/
	double fQR;			/*Activation energy for structural relaxation*/
	double ftauR0;		/*reference retardation time for structural relaxation*/
	double ftauRL, ftauRH; /*high and low limits of retardation time*/ 	
	
	double fTg;			/* glass transition temperature*/
	double fC1, fC2;	/*WLF constants*/
	double ftaug;		/*retardation time at Tg*/
	
	/*high temp stress response*/
//	double fmuN;		/*network stiffness*/
//	double flambdaL;	/*locking stretch*/
//	double fkappa;			/*bulk modulus*/
//	double fmueq;
	
	/*low temp viscous flow*/
//	double fmuneq;		/*nonequilibrium stiffness*/
	double fetaS0;		/*reference shear viscosity*/
//	double fetaSL, fetaSH;		/*high and low limits of etaSR*/
	double fQS;			/*activation energy for the viscoplastic flow*/
	double fsy0;		/*initial yield strength*/;
	double fsinf;		/*steady state limit of the yield strength*/
	double fh;			/*hardening modulus*/
		
	/*accessors for internal variables*/
	double* fsy;
	double* fsy_n;
	double* fdelneq;
	double* fdelneq_n;
	
	/*Residual*/
	dArrayT fRes;
	dArrayT fDelta;
	dMatrixT fGAB;
	dMatrixT fDAB;
	dMatrixT fDABbar;
	dMatrixT fMat;
};
}
#endif /* _SMP_simple_ */
