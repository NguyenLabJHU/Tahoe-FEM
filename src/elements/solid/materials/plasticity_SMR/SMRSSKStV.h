/* $Id: SMRSSKStV.h,v 1.1 2006-07-27 13:20:08 kyonten Exp $ */
/* created: Majid T. Manzari (04/16/2003) */
#ifndef _SMR_SS_KSTV_H_
#define _SMR_SS_KSTV_H_

/* base classes */
#include "SSIsotropicMatT.h"
#include "HookeanMatT.h"

namespace Tahoe {

/* forward declarations */
class SMRSSNLHardT;

class SMRSSKStV: public SSIsotropicMatT, public HookeanMatT
{
  public:

	/* constructor */
	SMRSSKStV(void);
	
	/* destructor */
	~SMRSSKStV(void);

	/* form of tangent matrix (symmetric by default) */
	virtual GlobalT::SystemTypeT TangentType(void) const;
	
	/** access strains from previous time step */
	//virtual bool Need_Strain_last(void) const { return true; };
	
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

	/** \name spatial description */
	/*@{*/
	/** spatial tangent modulus */
	virtual const dMatrixT& c_ijkl(void);
	virtual const dMatrixT& c_perfplas_ijkl(void);

	/** Cauchy stress */
	virtual const dSymMatrixT& s_ij(void);

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
 
  private:
  
  	/** pressure sensitive plasticity with nonlinear hardening and localization*/
	SMRSSNLHardT* fSMR;
  
  	/* return values */
  	dSymMatrixT	fStress;
  	dMatrixT	fModulus, fModulusCe;
    dMatrixT    fModulusPerfPlas;

};

} // namespace Tahoe 
#endif /* _SMR_SS_KSTV_H_ */
