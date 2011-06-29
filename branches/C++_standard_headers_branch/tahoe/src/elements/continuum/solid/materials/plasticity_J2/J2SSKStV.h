/* $Id: J2SSKStV.h,v 1.12 2005-07-26 16:31:24 paklein Exp $ */
/* created: paklein (06/18/1997) */
#ifndef _J2_SS_KSTV_H_
#define _J2_SS_KSTV_H_

/* base classes */
#include "SSSolidMatT.h"
#include "IsotropicT.h"
#include "HookeanMatT.h"
#include "J2SSC0HardeningT.h"

namespace Tahoe {

/** small strain J2 plastic material */
class J2SSKStV: public SSSolidMatT,
				public IsotropicT,
				public HookeanMatT,
				public J2SSC0HardeningT
{
public:

	/** constructor */
	J2SSKStV(void);

	/** model has history variables */
	virtual bool HasHistory(void) const { return true; };

	/* update internal variables */
	virtual void UpdateHistory(void);

	/* reset internal variables to last converged solution */
	virtual void ResetHistory(void);
	
	/** \name spatial description */
	/*@{*/
	/** spatial tangent modulus */
	virtual const dMatrixT& c_ijkl(void);

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
	virtual int NumOutputVariables(void) const;
	virtual void OutputLabels(ArrayT<StringT>& labels) const;
	virtual void ComputeOutput(dArrayT& output);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
	
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/* set modulus */
	virtual void SetModulus(dMatrixT& modulus);

        /*internal variable accessors*/
        virtual const iArrayT& InternalDOF(void) const;
        virtual const dArrayT& InternalStrainVars(void);
        virtual const dArrayT& InternalStressVars(void);

private:

	/* return values */
	dSymMatrixT	fStress;
	dMatrixT	fModulus;

	/** number of elastic iterations */
	int fElasticIterations;
};

} // namespace Tahoe 
#endif /* _J2_SS_KSTV_H_ */
