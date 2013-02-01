/* $Id: SMP_multisolvent.h,v 1.1 2013-02-01 17:00:11 tahoe.xiaorui Exp $ */
/* created: RX(08/01/2012) */
#ifndef _SMP_multisolvent_
#define _SMP_multisolvent_

/* base class */
#include "LAdMatrixT.h"
#include "RGSplitT2.h"
#include "InvLangevin.h"
#include "Gamma.h"

namespace Tahoe {

class SMP_multisolvent: public RGSplitT2
{
   public:

	/* constructor/destructor */
	SMP_multisolvent(void);

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

	/*These are the main functions*/
	/*stress response*/
	/*returns the free energy density the current deformation and internal state variables*/
	virtual double StrainEnergyDensity(void); 
	/*returns the spatial tangent modulus in reduced matrix form*/
	virtual const dMatrixT& c_ijkl(void);
	/*returns the symmetric Cauchy stress tensor*/
	virtual const dSymMatrixT& s_ij(void);
	/*returns inverse(F_T)*/
	virtual const dMatrixT& ThermalDeformation_Inverse(void);
	
	
   private:
	/*fictive temperature*/
	/*returns the fictive temperature Tf given the departure from equilibrium delta_bar_neq.*/
	double FictiveTemperature(const double deltaneq);
	/* return the time */
	virtual double Get_temp(double& timet);
	/*relaxation spectrum*/
	/*returns  tau_R/tauR0, the normalized strutural relaxation time*/
	double StructuralRelaxationFunc(const double Temperature, const double deltaneq, double theta);
	/*returns  tau_S/tauS0, normalized stress relaxation time*/
	double StressRelaxationFunc(const double Temperature, const double deltaneq, const double smag, const double sy, double theta);	
	/*returns  tau_Y/tauY0, normalized  relaxation time for evolution of yield strength*/
	double YieldRelaxationFunc(const double Temperature, const double deltaneq, const double smag, const double sy, double theta);	
	
	/*returns the relaxation spectrum for x=tau/taukWW*/
//	double StretchedExponentialSpectrum(const double x, const double tauKWW, const double betakWW);
	
	/*calculate discrete spectrum*/
//	void Compute_discrete_spectrum(const double beta, const double tauWWW, const double taumin, const double taumax, const int num,
//		dArrayT& times, dArrayT& spectrum, StringT filename);
	/*Numerical integration of the evolution equations for the eigenvalues of the elastic deformation tensor be, and delta_bar_neq*/
	void Compute_delneq(const dArrayT& delneq_n, dArrayT& delneq);
	void Compute_le(const ArrayT<dSymMatrixT>& C_vn, ArrayT<dSymMatrixT>& C_v, const double& sy_n, double& sy);
	void Compute_Kneq(dMatrixT& Modulus);

   protected:
   /*Gamma Function*/
   Gamma fGamma;
   
	/*Reference Temperature*/
	double fT0; /*Initial temperature*/
	double fTg;			/* glass transition temperature*/
	//double fC1, fC2;	/*WLF constants*/
	double fT2;         /*Kauzmann temperature 8*/
	double fA1, fb1, fbeta; /*fA is the energy barieer, fbb is the hyperbolic coefficent of the heat capacity and fbeta is the parameter related with the solvent */

	
	/*thermal expansion*/
	double falphar;		/*the rubbery, high temperature CTE*/
	double falphag;		/*the glassy, low temperature CTE*/

	/*thermal expansion*/
	double fmur;		/*the rubbery, high temperature shear modulus*/
	double fmug;		/*the glassy, low temperature shear modulus*/
	
	/*low temp viscous flow*/
	double fQS;			/*activation energy for the viscoplastic flow*/
	double fsy0;		/*initial yield strength*/;
	double fsinf;		/*steady state limit of the yield strength*/
	double ftauY0;			/*hardening modulus*/
	
	/*discrete structural relaxation spectrum*/
	int fNumR;	
	dArrayT ftimesR, fdalpha;  
	StringT fInputR;  
	
	int fNumT;	
	dArrayT ftimesT, ftempT; 
	StringT fInputRR; 

	/*discrete stress relaxation spectrum*/
	int fNumS;
	dArrayT ftimesS, fdmu;
	StringT fInputS;  
	
//	int fNumT;	
//	dArrayT ftimesT, ftempT;  
	StringT fInputT; 
	
	/*accessors for internal variables*/
	double* fsy;
	double* fsy_n;
	dArrayT fdelneq;
	dArrayT fdelneq_n;
	double* ftime;
	double* ftime_n;
	
	double ftheta0;
	
	dArrayT	fl_tr;
	dArrayT	fle;
	dArrayT fsig;
	
	LAdMatrixT fKdel;
	dArrayT fRdel;
	
	LAdMatrixT fKAB;
	LAdMatrixT fKAB2;
	dArrayT fRes;
	
	dArrayT fGA0;
	dArrayT fGA1;
	dArrayT fGA2;
	dMatrixT fDAB;
	dMatrixT fMat;
	dMatrixT fCalg2;
};
}
#endif /* _SMP_multisolvent_ */
