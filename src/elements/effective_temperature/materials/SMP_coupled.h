/* $Id: SMP_coupled.h,v 1.1 2015-09-20 03:48:45 tahoe.vickynguyen Exp $ */
/* created: RX (11/01/2013) / */
#ifndef _SMP_coupled_
#define _SMP_coupled_

/* base class */
//#include "RGSplitT2.h"
#include "FSThermoMechMatT.h"
#include "InvLangevin.h"
#include "FSSolidMatT.h"

namespace Tahoe {

    class SMP_coupled: public FSThermoMechMatT //,public RGSplitT2
{
   public:

	enum EnergyType {kMooneyRivlin=0,
					kLangevin=1};
					 
	enum ViscType {kNone = -1, 
		       kSimple=0, 
					kPower=1,
					kBergStromBoyce=2}; 
  
	/* constructor/destructor */
	SMP_coupled(void);
    virtual double Pressure(void) const {
		return 0.0;
	};

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
    /*The physical meaning is dS/dT */
    virtual const dSymMatrixT& d_ij(void);
    /*the d(themrmal body force)/dT (d(rho*deltac*dotTe-Wp)/dC) */
    virtual const dSymMatrixT& h_ij(void);
    /*the spatial heat flux */
    virtual  const dArrayT& q_i(void);
	
	virtual double RetardationTime(const double Temperature,const double Te);
	/*returns  eta_S, the shear viscosity*/
	virtual double ShearViscosity(const double Temperature,const double Te, const double smag);
   
    /*b1 is the d(themrmal body force)/dT as (d(rho*deltac*dotTf-Wp)/dT) */
   virtual double b1(void);
    /*this is the residual for the thermal part represented as:rho*cg*dotT+rho*deltac*dotTe-Wp */
    virtual double heatres(void);
    /* Also incorporate the anistropic heat conduction law */
    
    virtual const dMatrixT& ThermalDeformation_Inverse(void);
    
    virtual const dArrayT& Effective_Temperature(void);
    virtual const dArrayT& Effective_Temperature_last(void);

	protected:
	virtual void Initialize(void);
    
    const dMatrixT& c_ijkleq(void);
    const dMatrixT& c_ijklneq(void);
    
    const dSymMatrixT& s_ijeq(void);
    const dSymMatrixT& s_ijneq(void);
    const FSThermoMechSupportT* fFSThermoMechSupport;
 
    
    
private:
	/*Computes the algorithmic tangent for the nonequilibrium stress response*/
	virtual void Compute_Calg(const dArrayT& tau_dev, const dSymMatrixT& dtau_dev, const double& tau_m,
                              const double& dtau_m, dMatrixT& Calg, dArrayT& dalg, const int type);
	/*Numerical integration of the evolution equations for the eigenvalues of the elastic deformation tensor be, and delta_bar_neq*/
	virtual void ComputeEigs_e(const dArrayT& eigenstretch, dArrayT& eigenstretch_e,
                               dArrayT& eigenstress, dSymMatrixT& eigenmodulus, const int process_num);
/*TDN: Needs comments*/
    double a1(const double Temperature, const double Te, const double smag);
    double a2(const double Temperature, const double Te, const double smag);
    double b2(const double Temperature, const double Te, const double smag);
   protected:

	/*Reference Temperature*/
	 /*Temperature at which rubbery elastic properties are measured, not to be confused with the initial temperature at time = 0.*/
	double fT0;
    
	/*thermal expansion*/
	double falphar;		/*the rubbery, high temperature CTE*/
	double falphag;		/*the glassy, low temperature CTE*/
	double ftauR0;		/*reference retardation time for structural relaxation*/
	double ftauRL, ftauRH; /*high and low limits of retardation time*/
	
    double fdeltac;
	double fTg;			/* glass transition temperature*/
	double fC1, fC2;	/*WLF constants*/
	double ftaug;		/*retardation time at Tg*/
    double fA;
	
	/*low temp viscous flow*/
	double fetaS0;		/*reference shear viscosity*/
	double fetaSL, fetaSH;		/*high and low limits of etaSR*/
	double fQS;			/*activation energy for the viscoplastic flow*/
	double fsy0;		/* yield strength*/;

    dArrayT fTe;
    dArrayT fTe_n;
	
	/*Residual*/
	dArrayT fRes;
	dArrayT fDelta;
    dArrayT fdalg;
	dMatrixT fGAB;
	dMatrixT fDAB;
	dMatrixT fDABbar;
    dArrayT fDABbar2;
	dMatrixT fMat;
    dMatrixT fMat2;
    dSymMatrixT Hij;
    dSymMatrixT fStress3D_EQ;
    dSymMatrixT fStress3D_NEQ;
    dMatrixT fneqstress;
    dSymMatrixT fStress_EQ;
    dSymMatrixT fStress_NEQ;
    dMatrixT    fModulus3D_EQ;
    dMatrixT    fModulus3D_NEQ;
    dMatrixT    fModulus_EQ;
    dMatrixT    fModulus_NEQ;
	
	/* spectral operations */
//    SpectralDecompT fSpectralDecompSpat;
};
}
#endif /* _SMP_simple_ */
