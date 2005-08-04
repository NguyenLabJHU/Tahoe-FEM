/* $Id: GRAD_MRSSKStV.h,v 1.9 2005-08-04 21:49:33 kyonten Exp $ */
/* created: Karma Yonten (03/04/2004)                   
   Gradient Enhanced MR Model
*/
#ifndef _GRAD_MR_SS_KSTV_H_
#define _GRAD_MR_SS_KSTV_H_

#include "dSymMatrixT.h"

/* base classes */
#include "MFGPSSSolidMatT.h"
#include "IsotropicT.h"
#include "HookeanMatT.h"

namespace Tahoe 
{

/* forward declarations */
class GRAD_MRSSNLHardT;

class GRAD_MRSSKStV: public MFGPSSSolidMatT, public IsotropicT, public HookeanMatT
{
  public:

	/* constructor */
	GRAD_MRSSKStV(void);
	
	/* destructor */
	~GRAD_MRSSKStV(void);
	
	/* form of tangent matrix (symmetric by default) */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/** model has history variables */
	virtual bool HasHistory(void) const { return true; };

	/* update internal variables */
	virtual void UpdateHistory(void);

	/* reset internal variables to last converged solution */
	virtual void ResetHistory(void);
	
	/** returns elastic strain (3D) */
	virtual const dSymMatrixT& ElasticStrain(
                const dSymMatrixT& totalstrain,
				const ElementCardT& element, int ip);
				
	/** returns Laplacian of elastic strain (3D) */
	virtual const dSymMatrixT& LapElasticStrain(
                const dSymMatrixT& lap_totalstrain, 
				const ElementCardT& element, int ip);

	/** \name spatial description */
	/*@{*/
	/** spatial tangent modulus */
	virtual const dMatrixT& c_ijkl(void);

	virtual const dMatrixT& c_perfplas_ijkl(void);
	
	/*@{*/
	virtual const dMatrixT& c_UU1_ijkl(void);
	virtual const dMatrixT& c_UU2_ijkl(void);
	virtual const dMatrixT& c_ULam1_ij(void);
	virtual const dMatrixT& c_ULam2_ij(void);
	virtual const dMatrixT& c_LamU1_ij(void);
	virtual const dMatrixT& c_LamU2_ij(void);
	virtual const dMatrixT& c_LamLam1(void);
	virtual const dMatrixT& c_LamLam2(void);
	/*@{*/

	/** Cauchy stress */
	virtual const dSymMatrixT& s_ij(void);
	
	/** yield function */
	// return the yield function to form the RHS of the 
	// consistency equation
	virtual const double& YieldF(void);

	/** return the pressure associated with the last call to 
	 * SolidMaterialT::s_ij. See SolidMaterialT::Pressure
	 * for more information. */
	virtual double Pressure(void) const { return fStress.Trace()/3.0; };
	/*@}*/

	/* returns the strain energy density for the specified strain */
	virtual double StrainEnergyDensity(void);
	
	/* returns the number of variables computed for nodal extrapolation
	 * during for element output, ie. internal variables */
	virtual int  NumOutputVariables(void) const;
	virtual void OutputLabels(ArrayT<StringT>& labels) const;
	virtual void ComputeOutput(dArrayT& output);
	
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
	/*@}*/

protected:

	/* set modulus */
	virtual void SetModulus(dMatrixT& modulus); 
	int loccheck;
	
	/*dSymMatrixT	fStrain_IP, fLapStrain_IP;
    dArrayT fLambdaPM_IP, fLapLambdaPM_IP;
    ElementCardT curr_element;
    int curr_ip, num_ip;*/
 
private:
  
	/** pressure sensitive plasticity with nonlinear hardening and localization*/
	GRAD_MRSSNLHardT* fGRAD_MR;
  
  	/* return values */
  	dSymMatrixT	fStress;
  	dMatrixT	fModulus, fModulusCe;
    dMatrixT    fModulusPerfPlas;
    dMatrixT    fModulusUU1, fModulusUU2;
    dMatrixT    fModulusULam1, fModulusULam2;
    dMatrixT    fModulusLamU1, fModulusLamU2;
    dMatrixT    fModulusLamLam1, fModulusLamLam2;
    double      fYieldFunction; 

};

} // namespace Tahoe 
#endif /* _GRAD_MR_SS_KSTV_H_ */
