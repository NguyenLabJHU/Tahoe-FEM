/* $Id: DPSSKStVLoc.h,v 1.4 2004-07-15 08:28:56 paklein Exp $ */
/* created: myip (06/01/1999) */
#ifndef _DP_SS_KSTV_LOC_H_
#define _DP_SS_KSTV_LOC_H_

/* base classes */
#include "SSSolidMatT.h"
#include "IsotropicT.h"
#include "HookeanMatT.h"
#include "DPSSLinHardLocT.h"

namespace Tahoe {

class DPSSKStVLoc: public SSSolidMatT,
				public IsotropicT,
				public HookeanMatT,
				public DPSSLinHardLocT
{
public:

	/* constructor */
	DPSSKStVLoc(ifstreamT& in, const SSMatSupportT& support);

	/* initialization */
	virtual void Initialize(void);

	/* form of tangent matrix (symmetric by default) */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/* update internal variables */
	virtual void UpdateHistory(void);

	/* reset internal variables to last converged solution */
	virtual void ResetHistory(void);

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

	/*
	* Test for localization using "current" values for Cauchy
	* stress and the spatial tangent moduli. Returns 1 if the
	* determinant of the acoustic tensor is negative and returns
	* the normal for which the determinant is minimum. Returns 0
	* of the determinant is positive.
	*/
	// not used; see ComputeOutput and DetCheckT
	int IsLocalized(dArrayT& normal);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
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
  
	/* return values */
	dSymMatrixT	fStress;
	dMatrixT fModulus, fModulusCe;
	dMatrixT fModulusPerfPlas;

};

} // namespace Tahoe 
#endif /* _DP_SS_KSTV_LOC_H_ */
