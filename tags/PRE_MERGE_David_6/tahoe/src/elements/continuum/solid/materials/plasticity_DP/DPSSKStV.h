/* $Id: DPSSKStV.h,v 1.16 2006-06-30 15:13:22 regueiro Exp $ */
/* created: myip (06/01/1999) */
#ifndef _DP_SS_KSTV_H_
#define _DP_SS_KSTV_H_

/* base classes */
#include "SSIsotropicMatT.h"
#include "HookeanMatT.h"

namespace Tahoe {

/* forward declarations */
class DPSSLinHardT;

class DPSSKStV: public SSIsotropicMatT,
				public HookeanMatT
{
  public:

	/** constructor */
	DPSSKStV(void);

	/** constructor */
	~DPSSKStV(void);


	/* required parameter flags */


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

	/** \name spatial description */
	/*@{*/
	/** spatial tangent modulus */
	virtual const dMatrixT& c_ijkl(void);
	
	/* elastic modulus */
	virtual const dMatrixT& ce_ijkl(void);

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

    virtual bool IsLocalized(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipdirs, 
							AutoArrayT <double> &detAs, AutoArrayT <double> &dissipations_fact);
	virtual bool IsLocalized(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipdirs, double &detA);
	virtual bool IsLocalized(AutoArrayT <dArrayT> &normals, AutoArrayT <dArrayT> &slipdirs);

protected:

	/* set modulus */
 	virtual void SetModulus(dMatrixT& modulus); 
 
private:

	/** Drucker-Prager plasticity with linear hardening */
	DPSSLinHardT* fDP;
  
  	/* return values */
  	dSymMatrixT	fStress;
  	dMatrixT	fModulus;

};

} // namespace Tahoe 
#endif /* _DP_SS_KSTV_H_ */
