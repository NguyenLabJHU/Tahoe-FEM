/* $Id: SMP_simple.h,v 1.6 2009-06-23 08:48:57 thao Exp $ */
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

	/*Bookkeeping: initialize, update/reset internal variables and the beginning and end of each time step */
	virtual void PointInitialize(void);              
	virtual void UpdateHistory(void); // element at a time
	virtual void ResetHistory(void);  // element at a time


	/** \name Input/Output: 
		Defines the model parameters
		Reads in model parameters input file
		Overloads ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;
	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*do nothing*/
	virtual void InitStep(void);

	/*define outputs*/
	/*compute output variables*/ 
	virtual int NumOutputVariables() const; 
	virtual void OutputLabels(ArrayT<StringT>& labels) const; 
	virtual void ComputeOutput(dArrayT& output);

	/*These are the main functions*/
	/*stress response*/
	/*returns the free energy density the current deformation and internal state variables*/
	virtual double StrainEnergyDensity(void); 
	/*returns the spatial tangent modulus in reduced matrix form*/
	virtual const dMatrixT& c_ijkl(void);
	/*returns the symmetric Cauchy stress tensor*/
	virtual const dSymMatrixT& s_ij(void);
	
	/*fictive temperature*/
	/*returns the fictive temperature Tf given the departure from equilibrium delta_bar_neq.*/
	virtual double FictiveTemperature(const double deltaneq);
	
	/*viscosity*/
	/*returns  tau_R, the strutural relaxation time*/
	virtual double RetardationTime(const double Temperature, const double deltaneq);
	/*returns  eta_S, the shear viscosity*/
	virtual double ShearViscosity(const double Temperature, const double deltaneq, const double smag, const double sy);	
	
	/**compute thermal strains*/
	/*returns inverse(F_T)*/
	virtual const dMatrixT& ThermalDeformation_Inverse(void);

	protected:
	virtual void Initialize(void);

   private:
	/* set inverse of thermal transformation - return true if active */
//	virtual bool SetInverseThermalTransformation(dMatrixT& F_trans_inv);  			
	/*Computes the algorithmic tangent for the nonequilibrium stress response*/
	virtual void Compute_Calg(const dArrayT& tau_dev, const dSymMatrixT& dtau_dev, const double& tau_m, 
						const double& dtau_m, dMatrixT& Calg, const int type);
	/*Numerical integration of the evolution equations for the eigenvalues of the elastic deformation tensor be, and delta_bar_neq*/
	virtual void ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e, 
	                   dArrayT& eigenstress, dSymMatrixT& eigenmodulus, const int process_num);
    
   protected:
	/*Reference Temperature*/
	double fT0; /*Temperature at which rubbery elastic properties are measured, not to be confused with the initial temperature at time = 0.*/
	
	/*thermal expansion*/
	double falphar;		/*the rubbery, high temperature CTE*/
	double falphag;		/*the glassy, low temperature CTE*/
	double fT2;
	double fQR;			/*Activation energy for structural relaxation*/
	double ftauR0;		/*reference retardation time for structural relaxation*/
	double ftauRL, ftauRH; /*high and low limits of retardation time*/ 	
	
	double fTg;			/* glass transition temperature*/
	double fC1, fC2;	/*WLF constants*/
	double ftaug;		/*retardation time at Tg*/
	
	/*low temp viscous flow*/
//	double fmuneq;		/*nonequilibrium stiffness*/
	double fetaS0;		/*reference shear viscosity*/
	double fetaSL, fetaSH;		/*high and low limits of etaSR*/
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
